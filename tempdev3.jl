
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
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);; 

# To Run:

d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
#fid = 3490795;
#location = 1;

locs_d = Dict( 1391330 => 90, 3490795 => 1) # dict of locations.

st_dict = Dict{Int64,NTuple{9,Int64}}()
StateRecord(locs_d, dyn, st_dict)   # dict of restricted states. 


d1[dyn.all[1].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[1].fid][StateKey(dyn.all[1], dyn.all[1].level)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 2)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 3)] = 0.0

# FIXME - finish writing the test for this.  


ExactChoice(d1, d2, locs_d, st_dict, 3490795, 1, p1, p2, dyn, false)

d1[dyn.all[1].fid]

"""
function ExactChoice(temp::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                     stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, 
                     nbs::Dict{Int64, Int64}, # this should be a {Fid, Loc} dict  
                     st_recs::Dict{Int64,NTuple{9,Int64}}, # these should contain restricted states.  
                     fid::Int64, 
                     location::Int64,
                     p1::patientcount,
                     p2::patientcount,
                     D::DynState;    
                     messages::Bool = true, 
                     β::Float64 = 0.95,
                     ϕ13::Float64 = 0.0, # FIXME - substitute correct values and scale.
                     ϕ12::Float64 = 0.0,
                     ϕ1EX::Float64 = 0.0,
                     ϕ23::Float64 = 0.0,
                     ϕ21::Float64 = 0.0,
                     ϕ2EX::Float64 = 0.0,
                     ϕ31::Float64 = 0.0,
                     ϕ32::Float64 = 0.0,
                     ϕ3EX::Float64 = 0.0)   
    cps::Dict{Int64,Array{Float64,1}} = ContProbs(fid, st_recs, stable)  
    nstates::Dict{NTuple{9,Int64},Float64} = TotalCombine(D, location, nbs, cps)
    CV1::Float64 = ContVal(nstates, fid, stable ,1)
    CV2::Float64 = ContVal(nstates, fid, stable ,2)
    CV3::Float64 = ContVal(nstates, fid, stable ,3)  
    # Adding a test for the incrementing NaN problem.  #NANFIX
    testfloat::Float64 = deepcopy(D.all[location].mk.m[1].putils[2, findfirst(D.all[location].mk.m[1].putils[1,:], D.all[location].fid)])
  # Update value at Level 1
    D.all[location].level = 1
    UpdateD(D.all[location])                                  # updates the utility for a new level 
    DSimNew( D.all[location].mk, fid, p1, p2) # NANFIX                 # Computes the demand for that level.
    temp[fid][NStateKey(st_recs[fid],1)] = maximum([ϕ1EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([β*(CV1),-ϕ12+β*(CV2),-ϕ13+β*(CV3)])])
    UtilDown(D.all[location])                                 # resets the utility and level
    PatientZero(p1, p2)                                       # overwrites the patientcount with zeros 
    # NANFIX 
    if !isapprox(testfloat, D.all[location].mk.m[1].putils[2, findfirst(D.all[location].mk.m[1].putils[1,:], D.all[location].fid)])
        error("First Update Failed.")
    end 
  # Update value at Level 2 (repeats steps above!)
    D.all[location].level = 2
    UpdateD(D.all[location]) # Updates deterministic part of utility.  
    DSimNew( D.all[location].mk, fid, p1, p2) # NANFIX # Computes the demand.  
    temp[fid][NStateKey(st_recs[fid],2)] = maximum([ϕ2EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ21+β*(CV1),β*(CV2),-ϕ23+β*(CV3)])])
    UtilDown(D.all[location])                                 # resets the utility and level
    PatientZero(p1, p2)
    # NANFIX 
    if !isapprox(testfloat, D.all[location].mk.m[1].putils[2, findfirst(D.all[location].mk.m[1].putils[1,:], D.all[location].fid)])
        error("Second Update Failed.")
    end 
  # Update value at Level 3
    D.all[location].level = 3
    UpdateD(D.all[location])
    DSimNew( D.all[location].mk, fid, p1, p2) # NANFIX # Computes the demand.
    temp[fid][NStateKey(st_recs[fid],3)] = maximum([ϕ3EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ31+β*(CV1),-ϕ32+β*(CV2),β*(CV3)])])
    UtilDown(D.all[location])                                 # resets the utility and level
    PatientZero(p1, p2)
    # NANFIX 
    if !isapprox(testfloat, D.all[location].mk.m[1].putils[2, findfirst(D.all[location].mk.m[1].putils[1,:], D.all[location].fid)])
        error("Third Update Failed.")
    end 
end 




#=


=#




