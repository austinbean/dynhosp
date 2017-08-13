"""
`InexactConvergence(D::DynState,
                    chunk::Array{Int64,1}, 
                    p1::patientcount,
                    p2::patientcount,
                    tracker::Dict{NTuple{10, Int64}, Int64}, 
                    outvals::Dict{Int64, Dict{NTuple{10, Int64},Float64}},
                    itlim::Int64)`
"""
function InexactConvergence(D::DynState,
                            chunk::Array{Int64,1}, 
                            p1::patientcount,
                            p2::patientcount,
                            tracker::Dict{NTuple{10, Int64}, Int64}, 
                            outvals::Dict{Int64, Dict{NTuple{10, Int64},Float64}},
                            itlim::Int64)
    approxvals::Dict{NTuple{10, Int64},Float64} =Dict{NTuple{10, Int64},Float64}()        # put results of approximation here.  
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
        for i = 1:itlim
            nextstate = DictRandomize(outvals[k], elts, wgts)                             # gets the next state. 
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

