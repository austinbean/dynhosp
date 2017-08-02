"""
`ExactVal(D::DynState, V::allvisits, itlim::Int64, chunk::Array{Int64,1}; itlim::Int64 = 1, outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }())`
Computes the exact solution for smaller markets.  1 - 5 firms at most.
Note that this can take an argument with existing values.  This can be restarted, in other words.  

entries to consideR: 1-15 are all.
Here number of neighbors and entry (dyn.all[x])
Duopoly:
1, 4, 5, 9

triopoly:
18, 20, 21, 24, 25, 26

4-opoly
13, 14, 15


# testing: 


dyn = CounterObjects(5);
ch = [1] # first element
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
out1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
ExactVal(dyn, ch, p1, p2; outvals = out1 )



ExactVal(dyn, [11], p1, p2; outvals = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }());

dyn = CounterObjects(5);
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
ch2 = [11] # larger market. 
out2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
ExactVal(dyn, ch2, p1, p2; itlim = 2, outvals = out2)

# The following debugger will compare utility values at each point during the sim, but it is very slow.  
Insert this line somewhere in ExactVal:   



"""
function ExactVal(D::DynState,
                  chunk::Array{Int64,1}, 
                  p1::patientcount,
                  p2::patientcount;
                  itlim::Int64 = 100,
                  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }())
  tempvals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  DictClean(tempvals)                                                                   # initialize to zero. 
  totest::Dict{Int64,Bool} = Dict{Int64,Bool}()                                         # will record convergence (Fid, Bool) Dict.
  all_locs::Dict{Int64,Int64} = Dict{Int64, Int64}()                                    # will record the locations of competitors (Fid, Loc in Dyn.all) Dict, NOT the firm itself.  
  # TODO - this is kind of dumb... why not record this as a mutable cns type?  Can be converted to Tuple for Dict.
  # Or return a tuple in GiveState.
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
    tempvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    CompsDict(neighbors, D, all_locs)                                                   # now this is Dict{Fid, Location}
    StateRecord(all_locs, D, st_dict)                                                   # returns the restricted state.
    for el2 in neighbors                                                                # adds keys for the neighbors to the temp dict. 
      outvals[D.all[el2].fid] = Dict{NTuple{10, Int64}, Float64}()
      tempvals[D.all[el2].fid] = Dict{NTuple{10, Int64},Float64}() 
      totest[D.all[el2].fid] = false                                                    # initialized to FALSE - not converged. 
      StateEnumerate(TupletoCNS(st_dict[D.all[el2].fid]), outvals[D.all[el2].fid]) 
      StateEnumerate(TupletoCNS(st_dict[D.all[el2].fid]), tempvals[D.all[el2].fid])
    end 
    if !haskey(outvals, D.all[el].fid)
      StateEnumerate(D.all[el].cns, outvals[D.all[el].fid])
      StateEnumerate(D.all[el].cns, tempvals[D.all[el].fid])                            # this does NOT need starting values.  
    end 
    if haskey(outvals, D.all[el].fid)
      DictCopy(tempvals, outvals, 1.0)                                                  # if there is an entry for the value, copy FROM outvals TO tempvals.  
    end                                
    totest[D.all[el].fid] = false                                                       # all facilities to do initially set to false.  
  end
  altstates = MakeStateBlock(nfds)                                                      # generates a list of states to try, e.g., entry, exit and levels for each possible competitor.  
  converge::Bool = true
  #d2 = CounterObjects(5) # TODO - REMOVE.

  while (converge)&(its<itlim)                                                          # if true keep going.  
    for k in keys(totest)                                                              
      if !totest[k]  
        for r in 1:size(altstates,1)                                                    # Chooses a configuration. NB: Iterating over rows is not a great idea 
          # Be careful about this next line... is it changing the state in the right places?  
          # TODO - GiveState is definitely not working.   
          st_dict[k] = GiveState( D, chunk, all_locs, altstates[r,:], D.all[all_locs[k]].cns) 
          MapCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:])
          #TODO - print the share implied by WTP.    
          ExactChoice(tempvals, outvals, all_locs, st_dict, k, all_locs[k], p1, p2, D; messages = false)
          #println("calling reset !")
          ResetCompState(D, all_locs, chunk, FindFids(D, chunk), altstates[r,:]) # set it back 
          #DynAudit(D, d2) # TODO - remove  
        end 
      end 
    end
    # Convergence Test - this modifies bools in totest.
    ExactConvergence(tempvals, outvals, totest, its; messages = false)   
    # FIXME - prints minimum difference every 1,000 iterations.  Should be removed eventually. 
    if its %10 == 0
      println("iteration: ", its)
      for k1 in keys(outvals)
        mins::Float64 = 1.0
        for k2 in keys(outvals[k1])
          if tempvals[k1][k2] > 0
            if abs(outvals[k1][k2] - tempvals[k1][k2]) < mins
              mins = abs(outvals[k1][k2] - tempvals[k1][k2])
            end
          end 
          #println(k1, " ", k2, " ", mins) # just print the minimum difference.  
        end 
      end 
    end 
    # Copy the values and clean up.
    # TODO - check the weighting system here... iteration numbers so far are really low.  
    #PureCopy(outvals, tempvals)
    DictCopy(outvals, tempvals, 1/(its+1))                                            # NB - weight placed on new vs. old values.  
    DictClean(tempvals)                                                               # sets up for rewriting.
    for ky1 in keys(totest)                                                           # this tests every facility every time, but that's ok. 
      converge = ConvTest(totest)                                                     # iterates over bools in totest returns product
    end   
    its += 1
    #println("converge? ", converge)
  end 
  # Return equilibrium values:
  return outvals
end




#=







=#
