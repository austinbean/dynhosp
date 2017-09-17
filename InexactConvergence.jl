"""
`InexactConvergence(D::DynState,
                    chunk::Array{Int64,1}, 
                    p1::patientcount,
                    p2::patientcount,
                    tracker::Dict{NTuple{10, Int64}, Int64}, 
                    outvals::Dict{Int64, Dict{NTuple{10, Int64},Float64}},
                    itlim::Int64)`

# generates the tracker to test with InexactConvergence
As input take the output of NewApprox.
function tgen(d1)
  tr = Dict{NTuple{10, Int64}, Int64}()
  for k1 in keys(d1)
    for k2 in keys(d1[k1])
      tr[k2] = 100
    end 
  end 
  return tr 
end 
dyn = CounterObjects(1);
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
ddd = NewApprox(dyn, [37], p1, p2; wlh = 0, wlm = 10, itlim = 100)

tr1 = tgen(ddd)
InexactConvergence(dyn, [37], p1, p2,  tr1, ddd, 1000; wh = 0, wm = 6)

Slowdown here comes from |states|×|itlim| iterations.

InexactConvergence(dyn, [37], p1, p2, tr1, ab, 10) 138.363867 seconds (349.86 M allocations: 24.028 GiB, 1.82% gc time)

InexactConvergence(dyn, [37], p1, p2, tr1, ab, 100) 1425.038687 seconds (3.50 G allocations: 240.145 GiB, 1.81% gc time)

"""
function InexactConvergence(D::DynState,
                            chunk::Array{Int64,1},
                            p1::patientcount,
                            p2::patientcount,
                            tracker::Dict{NTuple{10, Int64}, Int64}, 
                            outvals::Dict{Int64, Dict{NTuple{10, Int64},Float64}},
                            itlim::Int64; 
                            wh::Int64 = 100, 
                            wm::Int64 = 0)
    strt = now()
    wl = Dates.Millisecond(Dates.Hour(wh)) + Dates.Millisecond(Dates.Minute(wm)) - Dates.Millisecond(Dates.Minute(5))
    approxvals::Dict{NTuple{10, Int64},Float64} = Dict{NTuple{10, Int64},Float64}()       # put results of approximation here.  
    all_locs::Dict{Int64,Int64} = Dict{Int64, Int64}()                                    # will record the locations of competitors (Fid, Loc in Dyn.all) Dict, NOT the firm itself.  
    st_dict::Dict{Int64,NTuple{9,Int64}} = Dict{Int64,NTuple{9,Int64}}()                  # will record the states of all firms from the point of view of el.
    neighbors::Array{Int64,1} = Array{Int64,1}()                                          # will record the locations in D.all[] of competing firms AND the firm itself.
    nfds::Array{Int64,1} = Array{Int64,1}()                                               # records the fids of neighbors, as fids, not locations. 
    for el in chunk                                                                       # goal of this loop is to: set up the dictionaries containing values with entries for the fids.  
        FindComps(D, neighbors, D.all[el])                                                # these are addresses of competitors in D.all 
        NFids(D, nfds, D.all[el])                                                         # records the fids of neighbors only, as fids.  
        push!(neighbors, el)                                                              # add the location of the firm in chunk
        CompsDict(neighbors, D, all_locs)                                                 # now this is Dict{Fid, Location}
    end
    altstates = MakeStateBlock(nfds)                                                      # generates a list of states to try, e.g., entry, exit and levels for each possible competitor.  
    inv_costs = zeros(9)
    wgts = zeros(outvals[D.all[chunk[1]].fid].count)                                      # this will hold the weights.  
    elts = Array{NTuple{10,Int64}}(outvals[D.all[chunk[1]].fid].count)                    # this will hold the states
    dists = RecordDists(D, chunk, all_locs)                                               # for each neighbor, record the distance to the main firm.  Do this once.
    statehold = zeros(Int64,9)
    k::Int64 = D.all[chunk[1]].fid                                                        # track total number of visits made.  
    totv::Int64 = 0
    for k1 in keys(tracker)                                                               # these are states, complete with level.
        appr::Float64 = 0
        current = now()
        # Now this breaks out of the loop, but not the whole function!
        if (current-strt)>wl 
            println("Convergence check time exceeded!") 
            return 0.0 # maybe it's the "return" that's causing the break, not the "break" itself.  
            break 
        end 
        for i = 1:itlim
            nextstate = DictRandomize(outvals[k], elts, wgts)                             # gets the next state. 
            # TODO - scale the euler constant in App Continuation.  
            c::Float64 = AppContinuation(nextstate, outvals[k])           
            r::Int64 = MakeConfig(nextstate, dists, altstates, statehold)
            MapCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:])
            InvCosts(StateShorten(nextstate), false, inv_costs)
            v::Float64 = AppChoice(all_locs[k], k, p1, p2, D, inv_costs)                  # this has roughly the right scale.
            appr += c + v 
            ResetCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:])  
        end 
        if haskey(approxvals, k1)
            approxvals[k1] += appr/itlim
        else 
            approxvals[k1] = appr/itlim
        end 
        totv += tracker[k1]
    end 
    # compute test:
    test::Float64 = 0.0 
    for k1 in keys(tracker)
        test += tracker[k1]*(outvals[k][k1] - approxvals[k1])^2
    end 
    return test/totv
end 

