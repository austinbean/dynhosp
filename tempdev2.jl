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
    neighbors::Array{Int64,1} = FindComps(D.all[location], D) # find the competitors.  
    # FIXME 04/25/2017 - the call to ContProbs already fs up because its looking for something that isn't there. 
    # but it isn't looking for it in stable.  it's in nlocs. 
    println("Exact Choice called on: ", fid)
    nloc = Array{Int64,1}() # initialize empty.
    if nbs[D.all[location].fid] # this will be true when this is a "neighbor" only, not a real facility.  
      nloc = FindComps(D.all[location], D) 
    else
      nloc = Array{Int64,1}() # I guess this array will be empty when the firm is the neighbor.  
    end 
    # FIXME - there are two calls to this function.  I don't seem to use "recs"
    # Also this is getting neighbors of the neighbor, which I do not want. 
    # Here is the issue - this calls on neighbors.  It should not.   
    recs = StateRecord(neighbors, location, D) # this returns the dict of the state for the firm whose value is being computed from point of view of the main firm.
    rec = StateRecord(D.all[location].nfids, location, D)     # this is computing the state from the point of view of... the main fac. 
    println("From Exact Choice ") 
    println("keys of stable: ", keys(stable)) 
    println(" temp keys: ",keys(temp))
# 04/26/2017 FIXME - the extra fid is here.  This searches in ContProbs USING rec.  
    println("state record ", rec)   
    println("nlocs: ", nloc)
    for el in nloc 
      println("fid is: ", D.all[el].fid)
    end 
    println("the fid: ", fid)
    println("set of neighbors: ", neighbors)
    println("Exact Choice nbs: ", nbs)
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








