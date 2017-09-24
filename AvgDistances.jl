
mutable struct DR # this is... how many patients traveled what distances from a zip. 
    p::patientcount 
    d::Float64 
end 


"""
`AverageD()`
How much does the average patient travel in excess of the regular arrangement?

- For each zip in a market, compute the total number of miles traveled, plus the average (easy)
- Reallocate facilities, update utilities (easy)
- Now recompute distances - see how much they increase.
- What are the costs of that travel?


dyn = CounterObjects(1);
chunk = [1];
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

AverageD(dyn, chunk, p1, p1)
"""
function AverageD(D::DynState,
                  chunk::Array{Int64,1}, 
                  p1::patientcount,
                  p2::patientcount)
  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  totest::Dict{Int64,Bool} = Dict{Int64,Bool}()                                         # will record convergence (Fid, Bool) Dict.
  all_locs::Dict{Int64,Int64} = Dict{Int64, Int64}()                                    # will record the locations of competitors (Fid, Loc in Dyn.all) Dict, NOT the firm itself.  
  st_dict::Dict{Int64,NTuple{9,Int64}} = Dict{Int64,NTuple{9,Int64}}()                  # will record the states of all firms from the point of view of el.
  neighbors::Array{Int64,1} = Array{Int64,1}()                                          # will record the locations in D.all[] of competing firms AND the firm itself.
  nfds::Array{Int64,1} = Array{Int64,1}()                                               # records the fids of neighbors, as fids, not locations. 
  its::Int64 = 1                                                                        # records iterations, but will be dropped after debugging.
  for el in chunk                                                                       # goal of this loop is to: set up the dictionaries containing values with entries for the fids.  
    FindComps(D, neighbors, D.all[el])                                                  # these are addresses of competitors in D.all 
    NFids(D, nfds, D.all[el])                                                           # records the fids of neighbors only, as fids.  
    push!(neighbors, el)                                                                # add the location of the firm in chunk
    if !haskey(outvals, D.all[el].fid)                                                  # add an empty dict IF there isn't already an entry.
      outvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    end 
    CompsDict(neighbors, D, all_locs)                                                   # now this is Dict{Fid, Location}
    StateRecord(all_locs, D, st_dict)                                                   # returns the restricted state.
    for el2 in neighbors                                                                # adds keys for the neighbors to the temp dict. 
      outvals[D.all[el2].fid] = Dict{NTuple{10, Int64}, Float64}()
      totest[D.all[el2].fid] = false                                                    # initialized to FALSE - not converged. 
      StateEnumerate(TupletoCNS(st_dict[D.all[el2].fid]), outvals[D.all[el2].fid]) 
    end                                 
    totest[D.all[el].fid] = false                                                       # all facilities to do initially set to false. 
    push!(nfds, D.all[el].fid) 
  end
  # need to keep facility specific patient counts.  
  # For each hospital (fid), how many people traveled how far to the hospital (from each zip)
  medcounts = Dict{NTuple{10,Int64}, Dict{Int64,Array{DR,1} } }()
  privcounts = Dict{NTuple{10,Int64}, Dict{Int64, Array{DR,1}}}()
  totald = Dict{Int64,Float64}()
  pcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  mcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  temparr = zeros(2, 12)
  altstates = MakeStateBlock(nfds)                                                      # generates a list of states to try, e.g., entry, exit and levels for each possible competitor.  
  converge::Bool = true
  inv_costs = zeros(9)
  k = D.all[chunk[1]].fid  
  while (its<2)                                                                         # if true keep going.  
    # just do this for the key of the main firm.  
    # create the state - and do this for ALL levels. 
    # This is not right - we need to compute ALL demands for ALL facilities.  But that's quite doable.
    # Enumerate facilities, list all patients, run this demand est on each hosp,
    # at each zip, record distance, this multiplied by numbers gives total miles.  Then divide by number of people.
    # Can incorporate mortality here too.                                                        
    for r in 1:size(altstates,1)                                                      # Chooses a configuration. 
      st_dict[k] = GiveState( D, chunk, all_locs, altstates[r,:], D.all[all_locs[k]].cns) 
      MapCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:])
      original = D.all[chunk[1]].level                          # save the orginal level. 
      D.all[chunk[1]].level = 1
      UpdateD(D.all[all_locs[k]])                                  # updates the utility for a new level only for the main firm.
      k1 = NStateKey(st_dict[k],1) 
      medcounts[k1] =  Dict{Int64,Array{DR,1} }()
      privcounts[k1] =  Dict{Int64,Array{DR,1} }()
      for k in keys(totest) # this will go over all of the neighboring firms. 
        medcounts[k1][D.all[all_locs[k]].fid] =  Array{DR,1}()
        privcounts[k1][D.all[all_locs[k]].fid] =  Array{DR,1}()
        for i = 1:size(D.all[all_locs[k]].mk.m,1)
            # these compute the demands.
            DemComp(D.all[all_locs[k]].mk.m[i].putils, temparr, pcount, D.all[all_locs[k]].fid, PatExpByType(D.all[all_locs[k]].mk.m[i].pcounts, true))  # pcount is an empty, pre-allocated patientcount into which results are written.
            DemComp(D.all[all_locs[k]].mk.m[i].mutils, temparr, mcount, D.all[all_locs[k]].fid, PatExpByType(D.all[all_locs[k]].mk.m[i].mcounts, false)) 
            # copy the pcount and mcount 
            mc1, pc1 = CopyCount(pcount, mcount) 
            # now measure the distance.
            d1 = distance(D.all[all_locs[k]].lat, D.all[all_locs[k]].long, D.all[all_locs[k]].mk.m[i].lat, D.all[all_locs[k]].mk.m[i].long)
            # push the copied p and m 
            push!(medcounts[k1][D.all[all_locs[k]].fid], DR(mc1, d1))
            push!(privcounts[k1][D.all[all_locs[k]].fid], DR(pc1, d1))
            # reset patient counts 
            ResetP(pcount)
            ResetP(mcount) 
        end 
      end 
      UtilDown(D.all[all_locs[k]])
      PatientZero(pcount, mcount)
      ResetCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:]) # set it back 
      # Level 2
      k2 = NStateKey(st_dict[k],2)
      D.all[chunk[1]].level = 1
      UpdateD(D.all[all_locs[k]])
      medcounts[k2] =  Dict{Int64,Array{DR,1} }()
      privcounts[k2] =  Dict{Int64,Array{DR,1} }()
      for k in keys(totest) # this will go over all of the neighboring firms. 
        medcounts[k2][D.all[all_locs[k]].fid] =  Array{DR,1}()
        privcounts[k2][D.all[all_locs[k]].fid] =  Array{DR,1}()
        for i = 1:size(D.all[all_locs[k]].mk.m,1)
            # these compute the demands.
            DemComp(D.all[all_locs[k]].mk.m[i].putils, temparr, pcount, D.all[all_locs[k]].fid, PatExpByType(D.all[all_locs[k]].mk.m[i].pcounts, true))  # pcount is an empty, pre-allocated patientcount into which results are written.
            DemComp(D.all[all_locs[k]].mk.m[i].mutils, temparr, mcount, D.all[all_locs[k]].fid, PatExpByType(D.all[all_locs[k]].mk.m[i].mcounts, false)) 
            # copy the pcount and mcount 
            mc1, pc1 = CopyCount(pcount, mcount) 
            # now measure the distance.
            d1 = distance(D.all[all_locs[k]].lat, D.all[all_locs[k]].long, D.all[all_locs[k]].mk.m[i].lat, D.all[all_locs[k]].mk.m[i].long)
            # push the copied p and m 
            push!(medcounts[k2][D.all[all_locs[k]].fid], DR(mc1, d1))
            push!(privcounts[k2][D.all[all_locs[k]].fid], DR(pc1, d1))
            # reset patient counts 
            ResetP(pcount)
            ResetP(mcount) 
        end 
      end 
      UtilDown(D.all[all_locs[k]])
      PatientZero(pcount, mcount)
      # Level 3
      k3 = NStateKey(st_dict[k],3)
      D.all[chunk[1]].level = 1
      UpdateD(D.all[all_locs[k]])  
      medcounts[k3] =  Dict{Int64,Array{DR,1} }()
      privcounts[k3] =  Dict{Int64,Array{DR,1} }()
      for k in keys(totest) # this will go over all of the neighboring firms. 
        medcounts[k3][D.all[all_locs[k]].fid] =  Array{DR,1}()
        privcounts[k3][D.all[all_locs[k]].fid] =  Array{DR,1}()
        for i = 1:size(D.all[all_locs[k]].mk.m,1)
            # these compute the demands.
            DemComp(D.all[all_locs[k]].mk.m[i].putils, temparr, pcount, D.all[all_locs[k]].fid, PatExpByType(D.all[all_locs[k]].mk.m[i].pcounts, true))  # pcount is an empty, pre-allocated patientcount into which results are written.
            DemComp(D.all[all_locs[k]].mk.m[i].mutils, temparr, mcount, D.all[all_locs[k]].fid, PatExpByType(D.all[all_locs[k]].mk.m[i].mcounts, false)) 
            # copy the pcount and mcount 
            mc1, pc1 = CopyCount(pcount, mcount) 
            # now measure the distance.
            d1 = distance(D.all[all_locs[k]].lat, D.all[all_locs[k]].long, D.all[all_locs[k]].mk.m[i].lat, D.all[all_locs[k]].mk.m[i].long)
            # push the copied p and m 
            push!(medcounts[k3][D.all[all_locs[k]].fid], DR(mc1, d1))
            push!(privcounts[k3][D.all[all_locs[k]].fid], DR(pc1, d1))
            # reset patient counts 
            ResetP(pcount)
            ResetP(mcount) 
        end 
      end 
      UtilDown(D.all[all_locs[k]])
      PatientZero(pcount, mcount)
      D.all[chunk[1]].level = original
    end 
    its += 1
  end 
  return medcounts, privcounts
end

"""
`ResetP(pp::patientcount)`
Sets patientcount back to 0.
"""
function ResetP(pp::patientcount)
 pp.count385 = 0.0
 pp.count386 = 0.0
 pp.count387 = 0.0
 pp.count388 = 0.0
 pp.count389 = 0.0
 pp.count390 = 0.0
 pp.count391 = 0.0
end 

"""
`CopyCount(pc::patientcount, mc::patientcount)`
Copies the patientcounts.  Allocates new ones.  
"""
function CopyCount(pc::patientcount, mc::patientcount)
    p1::patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    p2::patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    p1.count385 += pc.count385 
    p1.count386 += pc.count386 
    p1.count387 += pc.count387 
    p1.count388 += pc.count388 
    p1.count389 += pc.count389 
    p1.count390 += pc.count390 
    p1.count391 += pc.count391 
    p2.count385 += mc.count385 
    p2.count386 += mc.count386 
    p2.count387 += mc.count387 
    p2.count388 += mc.count388 
    p2.count389 += mc.count389 
    p2.count390 += mc.count390 
    p2.count391 += mc.count391 
    return p1, p2 
end 


"""
`TakeAverage(mc::Dict{NTuple{10,Int64}, Dict{Int64,Array{DR,1} } }, pc::Dict{NTuple{10,Int64}, Dict{Int64, Array{DR,1}}})`
Compute the average distances traveled per market config.  
"""
function TakeAverage(mc::Dict{NTuple{10,Int64}, Dict{Int64,Array{DR,1} } }, pc::Dict{NTuple{10,Int64}, Dict{Int64, Array{DR,1}}}, f::Int64) 
  # outp... fid, nine states, one level, one count of patients, one average distance
  cols = 13
  rows = mc.count 
  outp = zeros(rows, cols)
  rc = 1 # row counter
  for k1 in keys(mc)       # the state config
    pats::Float64 = 0.0
    ds::Float64 = 0.0
    for k2 in keys(mc[k1]) # the other firms.
      for i = 1:size(mc[k1][k2],1)
        a1, b1 = DREX(mc[k1][k2][i])
        pats += a1
        ds += a1*b1
      end 
    end 
    println(ds)
    # now put it in the output in some way... 
    outp[rc,1] = f 
    for j = 1:length(k1)
      outp[rc, j+1] = k1[j]
    end 
    outp[rc,12] += pats 
    outp[rc,13] += (ds/pats)
    rc += 1
  end 
  return outp
end


"""
`DREX(d::DR)`
Returns sum of patients and distance.
"""
function DREX(d::DR)
  return sum(d.p), d.d 
end 




#=


=#
