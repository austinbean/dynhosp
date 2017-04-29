function ExactChoice(temp::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                     stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, 
                     nbs::Dict{Int64, Int64},
                     fid::Int64, 
                     location::Int64,
                     p1::patientcount,
                     p2::patientcount,
                     D::DynState; 
                     messages::Bool = true, 
                     β::Float64 = 0.95,
                     ϕ13::Float64 = 0.0, # substitute values and scale these!
                     ϕ12::Float64 = 0.0,
                     ϕ1EX::Float64 = 0.0,
                     ϕ23::Float64 = 0.0,
                     ϕ21::Float64 = 0.0,
                     ϕ2EX::Float64 = 0.0,
                     ϕ31::Float64 = 0.0,
                     ϕ32::Float64 = 0.0,
                     ϕ3EX::Float64 = 0.0)  
    # Also this is getting neighbors of the neighbor, which I do not want. 
    # use the Dict in nbs.    
#    recs = StateRecord(neighbors, location, D) # this returns the dict of the state for the firm whose value is being computed from point of view of the main firm.
    rec = StateRecord(D.all[location].nfids, location, D)     # this is computing the state from the point of view of... the main fac. 
    println("From Exact Choice ") 
    println("keys of stable: ", keys(stable)) 
    println(" temp keys: ",keys(temp))
    println("state record ", rec)    
    println("the fid: ", fid)
    println("set of neighbors: ", neighbors)
    println("Exact Choice nbs: ", nbs)


# FIXME - here there is nothing called nloc anymore.  

    cps::Dict{Int64,Array{Float64,1}} = ContProbs(rec, nloc, stable, D)  
    nstates::Dict{NTuple{9,Int64},Float64} = TotalCombine(D, location, D.all[location].nfids, cps)
    CV1::Float64 = ContVal(nstates, fid, stable ,1)
    CV2::Float64 = ContVal(nstates, fid, stable ,2)
    CV3::Float64 = ContVal(nstates, fid, stable ,3)   
    
    if messages 
      
      println("stable keys before", keys(stable))  
      println("CV's: ", CV1, " ", CV2, " ", CV3) 
      println("the max was ", temp[fid][StateKey(D.all[location],1)]) 
      println("the rev was ", PatientRev(D.all[location],p1,p2,10)) 
    end
  # Update value at Level 1
    D.all[location].level = 1
    UpdateD(D.all[location])                                  # updates the utility for a new level 
    DSimNew( D.all[location].mk, fid, p1, p2)                 # Computes the demand for that level. 
    temp[fid][StateKey(D.all[location],1)] = maximum([ϕ1EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([β*(CV1),-ϕ12+β*(CV2),-ϕ13+β*(CV3)])])
    D.all[location].level = D.all[location].actual            # resets the level 
    UtilDown(D.all[location])                                 # resets the utility
    PatientZero(p1, p2)                                       # overwrites the patientcount with zeros 
  # Update value at Level 2 (repeats steps above!)
    D.all[location].level = 2
    UpdateD(D.all[location]) # Updates deterministic part of utility.  
    DSimNew( D.all[location].mk, fid, p1, p2) # Computes the demand.   
    temp[fid][StateKey(D.all[location],2)] = maximum([ϕ2EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ21+β*(CV1),β*(CV2),-ϕ23+β*(CV3)])])
    D.all[location].level = D.all[location].actual
    UtilDown(D.all[location])
    PatientZero(p1, p2)
  # Update value at Level 3
    D.all[location].level = 3
    UpdateD(D.all[location])
    DSimNew( D.all[location].mk, fid, p1, p2) # Computes the demand.
    temp[fid][StateKey(D.all[location],3)] = maximum([ϕ3EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ31+β*(CV1),-ϕ32+β*(CV2),β*(CV3)])])
    D.all[location].level = D.all[location].actual
    UtilDown(D.all[location])
    PatientZero(p1, p2)
end 



#=

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
 Ok - the thing is that this must be done for Both facilities and their neighbors, but the notion of the state 
 for neighbors is different. This is important. 


##### TESTING ######
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients); 

# To Run:

d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
#fid = 3490795;
#location = 1;

d1[dyn.all[1].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[1].fid][StateKey(dyn.all[1], dyn.all[1].level)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 2)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 3)] = 0.0


ExactChoice(d1, d2, dyn.all[1].fid, 1, p1, p2,  dyn)
d1[dyn.all[1].fid]


EXTRA: 

                # FIXME - here is a problem.  Keys are being added to these dicts inconsistently.  Don't add all of them.  Why 
                # are these needed anyway? 
                # I want to *not* add these.  What will that break?   
                # if !haskey(stable, fid) # this should not be necessary when this is debugged.  
                #   stable[fid] = Dict{NTuple{10, Int64},  Float64 }()
                # end 
                # for el in keys(recs) # this adds a record for each of the (state,level) options.  They are put in the stable dict.  
                #   if !haskey(stable, el)
                #     stable[el] = Dict{NTuple{10,Int64}, Float64}()
                #   end
                #   if !haskey( stable[el], TAddLevel(recs[el], 1) )
                #     stable[el][TAddLevel(recs[el], 1)] = 0.5
                #   end 
                #   if !haskey( stable[el], TAddLevel(recs[el], 2) )
                #     stable[el][TAddLevel(recs[el], 2)] = 0.5
                #   end 
                #   if !haskey( stable[el], TAddLevel(recs[el], 3) )
                #     stable[el][TAddLevel(recs[el], 3)] = 0.5
                #   end 
                # end 

    # FIXME - I don't want to FindComps when this is a neighbor.
    # should I permanently take as arguments all of the neighbors I want to do?  
    # think about the 2-3 firm case... especially with non-overlapping sets of neighbors.  
"""


=#




