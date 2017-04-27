
function ExactVal(D::DynState,
                  chunk::Array{Int64,1},
                  p1::patientcount,
                  p2::patientcount;
                  messages::Bool = false,
                  itlim::Int64 = 100,
                  debug::Bool = true,
                  beta::Float64 = 0.95,
                  conv::Float64 = 0.0001)
  if messages println("from exactval") end
  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  tempvals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  totest::Dict{Int64,Bool} = Dict{Int64,Bool}()                                       # will record convergence 
  locs::Dict{Int64,Int64} = Dict{Int64, Int64}()                                      # will record the locations of competitors 
  nbs::Dict{Int64, Bool} = Dict{Int64, Bool}()
  its::Int64 = 0                                                                      # records iterations, but will be dropped after debugging.
  for el in chunk # goal of this loop is to: set up the dictionaries containing values with entries for the fids.  
    # this is doing it wrong?  Is el getting reassigned?
    println("first el: ", el)
    neighbors::Array{Int64,1} = FindComps(D.all[el], D)                               #  these are addresses of fids.
    println("Competitors are: ", neighbors) # but this is printed in the next line?
    println("main fid is: ", D.all[el].fid)
    for el2 in neighbors
      println("printing out fid of neighbors: ", D.all[el2].fid) # this should find 1.
    end 
    println("second el: ", el)
    outvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    tempvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    println("the el is: ", el)
    println("call to StateRecord: ", D.all[el].nfids) # this has too many elements in it? Or what? 
    # why does this return something with three keys?  WTF?  
    stdict = StateRecord(D.all[el].nfids, el, D)                                      # returns the restricted state.
    println("st dict: ", keys(stdict))
    for el2 in neighbors                                                              # adds keys for the neighbors to the temp dict. 
      println("the fid ", D.all[el2].fid)
      println("the keys of outvals: ", keys(outvals))
      println("the keys of tempvals: ", keys(tempvals))
      println("keys of stdict: ", keys(stdict))
      outvals[D.all[el2].fid] = Dict{NTuple{10, Int64}, Float64}()
      tempvals[D.all[el2].fid] = Dict{NTuple{10, Int64},Float64}() 
      totest[D.all[el2].fid] = true                                                   # don't test convergence of neighbors temporarily.  FIXME 
      println("err1")
      StateEnumerate(TupletoCNS(stdict[D.all[el2].fid]), outvals[D.all[el2].fid]) 
      StateEnumerate(TupletoCNS(stdict[D.all[el2].fid]), tempvals[D.all[el2].fid])
      println("err2")
      locs[D.all[el2].fid] = el2
      nbs[D.all[el2].fid] = true # is it a neighbor?  Yes.
    end 
    StateEnumerate(D.all[el].cns, outvals[D.all[el].fid])                             # TODO - starting values here.
    StateEnumerate(D.all[el].cns, tempvals[D.all[el].fid])                            # this does NOT need starting values.  
    totest[D.all[el].fid] = true                                                      # all facilities to do initially set to false.  
    locs[D.all[el].fid] = el                                                          # stores a fid,location value
    nbs[D.all[el].fid] = false   # is it a neighbor? No.  
  end
  # Updating process:
  converge = false
  while (!converge)&(its<itlim)                                                       # if true keep going.    
    converge = true                                                                   # reassign, to catch when it terminates.
    for k in keys(totest)                                                             # TODO - not updating the competitors.
      if totest[k]                                                                    # only run those for which true.
        if messages println("tempvals keys before: ", keys(tempvals)) end
        if messages println("outvals keys before: ", keys(outvals)) end
        if messages println("locs keys before: ", keys(locs)) end 
        ExactChoice(tempvals, outvals, nbs, k, locs[k], p1, p2, D; messages = true)  
        # if messages println("tempvals keys after: ", keys(tempvals)) end
        # if messages println("outvals keys after: ", keys(outvals)) end
        # if messages println("locs keys after: ", keys(locs)) end 
      end 
    end
    # Convergence Test:
    if messages println("near convergence") end 
    converge, totest = ExactConvergence(tempvals, outvals, totest; messages = true)    
    # Copy the values and clean up.
    DictCopy(outvals, tempvals)
    DictClean(tempvals)                                                               # sets up for rewriting.
    for ky1 in keys(totest)                                                           # this tests every facility every time, but that's ok. 
      converge = converge&totest[ky1]  
    end   
    its += 1
    println("iteration ", its)
  end 
  # Return equilibrium values:
  return outvals
end
