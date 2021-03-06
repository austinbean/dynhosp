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

dyn = pm.CounterObjects(1);
p1 = pm.patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = pm.patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
pm.NewApprox(dyn, [37], p1, p2; wlh = 0, wlm = 5, itlim = 200)


Try this approximation with 4916068 - > where is that one?  This one has eight neighbors.  Index is 37.  

# Testing remote call: 
remotecall_fetch(NewApprox, p, CounterObjects(1), [chs[ix]], patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0), patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0); wallh = 0, wallm = 2 )

# Problematic firms - these don't cancel in time.  Maybe this means they take more than 5 minutes to do 10 iterations?  Scarcely believable.
4916068, 3390720, 4536048, 4536337, 4536253, 293005, 293015, 296015, 296025, 293070, 296002, 293120
Likely problem: convergence test takes too long and does not check time.  




# Testing beginning:  
dyn = CounterObjects(1);
all_l = Dict{Int64,Int64}()
st_d = Dict{Int64, NTuple{9,Int64}}()
neb = Array{Int64,1}()
nfs = Array{Int64,1}()
outv = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
outv[dyn.all[52].fid] = Dict{NTuple{10, Int64}, Float64 }()

FindComps(dyn, neb, dyn.all[52])
NFids(dyn, nfs, dyn.all[52])
push!(neb, 52)
CompsDict(neb, dyn, all_l)
StateRecord(all_l, dyn, st_d)
StateEnumerate(dyn.all[52].cns, outv[dyn.all[52].fid])
altstates = MakeStateBlock(nfs)                                                        
inv_costs = zeros(9)
wgts = zeros(outv[dyn.all[52].fid].count)                                            
elts = Array{NTuple{10,Int64}}(outv[dyn.all[52].fid].count)                    
dists = RecordDists(dyn, [52], all_l)                                               
statehold = zeros(Int64,9)
k = dyn.all[52].fid

 @time NewApprox(dyn, [37], p1, p2; wlh = 0, wlm = 10, itlim = 1_0000)
 33.681478 seconds (94.45 M allocations: 6.483 GiB, 2.05% gc time)
Dict{Int64,Dict{NTuple{10,Int64},Float64}} with 1 entry:
  4916068 => Dict((0, 0, 0, 0, 2, 1, 1, 2, 0, 3)=>0.00343493,(0, 0, 0, 1, 2, 1,…


 @time NewApprox(dyn, [37], p1, p2; wlh = 0, wlm = 100, itlim = 1_000_000)
4198.062779 seconds (9.55 G allocations: 656.114 GiB, 1.82% gc time)
Dict{Int64,Dict{NTuple{10,Int64},Float64}} with 1 entry:
4916068 => Dict((0, 0, 0, 0, 2, 1, 1, 2, 0, 3)=>0.000258277,(0, 0, 0, 1, 2, 1, 2, 2, 0, 3)=>0.00015639…


"""
function NewApprox(D::DynState,
                   chunk::Array{Int64,1}, 
                   p1::patientcount,
                   p2::patientcount;
                   wlh::Int64 = 100, 
                   wlm::Int64 = 0,
                   itlim::Int64 = 100_000)
  strt = now()
  wl = Dates.Millisecond(Dates.Hour(wlh)) + Dates.Millisecond(Dates.Minute(wlm)) - Dates.Millisecond(Dates.Minute(10))
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
    # TODO - I think that this is causing a problem.  I am not sure that this state is the best.  
    StateRecord(all_locs, D, st_dict)                                                   # returns the restricted state. 
    StateEnumerate(D.all[el].cns, outvals[D.all[el].fid])
  end
  altstates = MakeStateBlock(nfds)                                                      # generates a list of states to try, e.g., entry, exit and levels for each possible competitor.  
  converge::Bool = true
  inv_costs = zeros(9)
  wgts = zeros(outvals[D.all[chunk[1]].fid].count)                                      # this will hold the weights.  
  elts = Array{NTuple{10,Int64}}(undef,outvals[D.all[chunk[1]].fid].count)              # this will hold the states
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
    ResetCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:])              # set it back 
    if its%10 == 0                                                                      # Check time every 100 iterations
      current = now()
      if (current-strt)>wl 
        println("Time exceeded!") 
        return outvals
        break 
      end 
    end 
    if its%1_000_000 == 0                                                               # Check convergence every 1,000,000 - takes 7 minutes for a large state.
      cv::Float64 = InexactConvergence(D, chunk, p1, p2, tracker, outvals, 100; wh = wlh, wm = wlm)
      if cv < 1e-6
        converge = false # stop.
      end 
      println("iteration: ", its, "  ", cv, " ", converge )
      ResetTracker(tracker) 
    end 
    its += 1
  end 
  # Return equilibrium values:
  return outvals
end
