"""
`NewApprox(D::DynState,
                  chunk::Array{Int64,1}, 
                  p1::patientcount,
                  p2::patientcount;
                  wallh::Int64 = 100, 
                  wallm::Int64 = 0,
                  itlim::Int64 = 100000,
                  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }())`

Second attempt at EP approximation.

dyn = CounterObjects(5);
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
ch2 = [11] # larger market. 
NewApprox(dyn, ch2, p1, p2; wallh = 0, wallm = 2, itlim = 200_000)


"""
function NewApprox(D::DynState,
                   chunk::Array{Int64,1}, 
                   p1::patientcount,
                   p2::patientcount;
                   wallh::Int64 = 100, 
                   wallm::Int64 = 0,
                   itlim::Int64 = 100_000)
  strt = now()
  wl = Dates.Millisecond(Dates.Hour(wallh)) + Dates.Millisecond(Dates.Minute(wallm)) - Dates.Millisecond(Dates.Minute(1))
  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  tracker::Dict{NTuple{10, Int64}, Int64} = Dict{NTuple{10, Int64}, Int64}()            # Dict to hold visited states
  all_locs::Dict{Int64,Int64} = Dict{Int64, Int64}()                                    # will record the locations of competitors (Fid, Loc in Dyn.all) Dict, NOT the firm itself.  
  st_dict::Dict{Int64,NTuple{9,Int64}} = Dict{Int64,NTuple{9,Int64}}()                  # will record the states of all firms from the point of view of el.
  neighbors::Array{Int64,1} = Array{Int64,1}()                                          # will record the locations in D.all[] of competing firms AND the firm itself.
  nfds::Array{Int64,1} = Array{Int64,1}()                                               # records the fids of neighbors, as fids, not locations. 
  its::Int64 = 1                                                                        # records iterations, but will be dropped after debugging.
  for el in chunk                                                                       # goal of this loop is to: set up the dictionaries containing values with entries for the fids.  
    FindComps(D, neighbors, D.all[el])                                                  # these are addresses of competitors in D.all 
    NFids(D, nfds, D.all[el])                                                           # records the fids of neighbors only, as fids.  
    push!(neighbors, el)                                                                # add the location of the firm in chunk
    outvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    CompsDict(neighbors, D, all_locs)                                                   # now this is Dict{Fid, Location}
    StateRecord(all_locs, D, st_dict)                                                   # returns the restricted state. 
    StateEnumerate(D.all[el].cns, outvals[D.all[el].fid])
  end
  altstates = MakeStateBlock(nfds)                                                      # generates a list of states to try, e.g., entry, exit and levels for each possible competitor.  
  converge::Bool = true
  inv_costs = zeros(9)
  wgts = zeros(outvals[D.all[chunk[1]].fid].count)                                      # this will hold the weights.  
  elts = Array{NTuple{10,Int64}}(outvals[D.all[chunk[1]].fid].count)                    # this will hold the states
  dists = RecordDists(D, chunk, all_locs)                                               # for each neighbor, record the distance to the main firm.  Do this once.
  statehold = zeros(Int64,9)
  k::Int64 = D.all[chunk[1]].fid 
  while (converge)&(its<itlim)                                                          # if true keep going.  
    nextstate = DictRandomize(outvals[k], elts, wgts)                                   # generate the next state.
    if haskey(tracker, nextstate)
      tracker[nextstate] += 1
    else 
      tracker[nextstate] = 1
    end 
    r = MakeConfig(nextstate, dists, altstates, statehold) 
    MapCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:])
    InvCosts(st_dict[k], false, inv_costs)
    outvals[k][nextstate] = (1-(1/(1+tracker[nextstate])))*(outvals[k][nextstate]) + (1/(1+tracker[nextstate]))*(AppChoice(all_locs[k], k, p1, p2, D, inv_costs)) # only map the state out to one element of outvals.
    #println("its: ", its, " min: ", OutMin(outvals, tracker, k))
    ResetCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:])              # set it back 
    if its %10_000 == 0
      cv::Float64 = InexactConvergence(D, chunk, p1, p2, tracker, outvals, 100)
      println(cv)
      if cv < 1e-6
        converge = false # stop.
      end 
      println(converge )
      ResetTracker(tracker)  
      println("iteration: ", its)
      current = now()
      if (current-strt)>wl 
        println("Time exceeded!") 
        return outvals
        break 
      end 
    end 
    its += 1
  end 
  # Return equilibrium values:
  return outvals
end
