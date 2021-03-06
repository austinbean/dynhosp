
"""
`ExactChoice(temp::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                     stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, 
                     fid::Int64, 
                     location::Int64,
                     p1::patientcount,
                     p2::patientcount,
                     competitors::Array{Int64,1},
                     D::DynState; )`
What action should the firm choose?
Takes two dictionaries, the DynState, computes the best action, returns the value of the action.
Needs to: 
- compute the demand in expectation at EACH possible level.
- compute the profit at EACH possible level.
- state will be recorded in the dyn record. 
- But the key thing is: return the VALUE of the state.   
NB - level won't change.  I can compute the value of being in all of these states depending on the level.



##### TESTING ######
TexasEq = pm.CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = pm.EntireState(Array{pm.Market,1}(), Dict{Int64, pm.Market}(), Dict{Int64, Int64}());
pm.CMakeIt(Tex, ProjectModule.fips);
pm.FillState(Tex, ProjectModule.alldists, 50);
patients = pm.NewPatients(Tex);

dyn = pm.DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);

# To Run:

d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = pm.patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = pm.patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
#fid = 3490795;
#location = 1;

locs_d = Dict( 1391330 => 90, 3490795 => 1) # dict of locations.

st_dict = Dict{Int64,NTuple{9,Int64}}()
pm.StateRecord(locs_d, dyn, st_dict)   # dict of restricted states. 


d1[dyn.all[1].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[1].fid][pm.StateKey(dyn.all[1], dyn.all[1].level)] = 0.0
d1[dyn.all[1].fid][pm.StateKey(dyn.all[1], 2)] = 0.0
d1[dyn.all[1].fid][pm.StateKey(dyn.all[1], 3)] = 0.0

# FIXME - finish writing the test for this.  


pm.ExactChoice(d1, d2, locs_d, st_dict, 3490795, 1, p1, p2, dyn, false)

d1[dyn.all[1].fid]

For debugging: 
UpdateUCheck(D.all[location]) 
can be used

"""
function ExactChoice(temp::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                     stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, 
                     nbs::Dict{Int64, Int64},              # this should be a {Fid, Loc} dict  
                     st_recs::Dict{Int64,NTuple{9,Int64}}, # TODO - change for NTuple{9,Int64} OR   
                     fid::Int64, 
                     location::Int64,
                     p1::patientcount,
                     p2::patientcount,
                     D::DynState,
                     inv_costs::Array{Float64};    
                     counter::Bool = false) 
    β::Float64 = 0.95
    ϕ13::Float64 = inv_costs[1] 
    ϕ12::Float64 = inv_costs[2]
    ϕ1EX::Float64 = inv_costs[3]
    ϕ23::Float64 = inv_costs[4]
    ϕ21::Float64 = inv_costs[5]
    ϕ2EX::Float64 = inv_costs[6]
    ϕ31::Float64 = inv_costs[7]
    ϕ32::Float64 = inv_costs[8]
    ϕ3EX::Float64 = inv_costs[9]   
    cps::Dict{Int64,Array{Float64,1}} = ContProbs(fid, st_recs, stable)  
    nstates::Dict{NTuple{9,Int64},Float64} = TotalCombine(D, location, nbs, cps)
    CV1::Float64 = ContVal(nstates, fid, stable ,1)
    CV2::Float64 = ContVal(nstates, fid, stable ,2)
    CV3::Float64 = ContVal(nstates, fid, stable ,3)  
  # Update value at Level 1
    original = D.all[location].level                          # save the orginal level. 
    D.all[location].level = 1
    UpdateD(D.all[location])                                  # updates the utility for a new level only for the main firm.
    DSimNew( D.all[location].mk, fid, p1, p2)                 # Computes the demand for that level.
    temp[fid][NStateKey(st_recs[fid],1)] = maximum([ϕ1EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([β*(CV1),-ϕ12+β*(CV2),-ϕ13+β*(CV3)])])
    UtilDown(D.all[location])                                 # resets the utility and level - just for the single firm.
    PatientZero(p1, p2)                                       # overwrites the patientcount with zeros 
  # Update value at Level 2  
    D.all[location].level = 2
    UpdateD(D.all[location])                                  # Updates deterministic part of utility for the main firm.
    DSimNew( D.all[location].mk, fid, p1, p2)                 # Computes the demand.  
    temp[fid][NStateKey(st_recs[fid],2)] = maximum([ϕ2EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ21+β*(CV1),β*(CV2),-ϕ23+β*(CV3)])])
    UtilDown(D.all[location])                                 # resets the utility and level for the main firm.  
    PatientZero(p1, p2)
  # Update value at Level 3
    D.all[location].level = 3
    UpdateD(D.all[location])
    DSimNew( D.all[location].mk, fid, p1, p2)                 # Computes the demand.
    temp[fid][NStateKey(st_recs[fid],3)] = maximum([ϕ3EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ31+β*(CV1),-ϕ32+β*(CV2),β*(CV3)])])
    UtilDown(D.all[location])                                 # resets the utility and level for the main firm.
    PatientZero(p1, p2)
    D.all[location].level = original                          # reset original level. 
end 




#=


=#




