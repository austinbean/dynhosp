function ExactVal(D::DynState,
                  chunk::Array{Int64,1}, 
                  p1::patientcount,
                  p2::patientcount;
                  messages::Bool = false,
                  itlim::Int64 = 100,
                  debug::Bool = true,
                  beta::Float64 = 0.95,
                  conv::Float64 = 0.0001)
  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  tempvals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  totest::Dict{Int64,Bool} = Dict{Int64,Bool}()                                       # will record convergence 
  locs::Dict{Int64,Int64} = Dict{Int64, Int64}()                                      # will record the locations of competitors 
  its::Int64 = 0                                                                      # records iterations, but will be dropped after debugging.
  for el in chunk # goal of this loop is to: set up the dictionaries containing values with entries for the fids.  
    # Now this will take ONE firm, that in "chunk", then use that to find the relevant neighbors, but NOT look for neighbors of those firms.
    # consider rewriting FindComps to return a Dict{Int64, Int64} = {Fid, Loc}.  To replace locs.  
    neighbors::Array{Int64,1} = FindComps(D, D.all[el])                             #  these are addresses of fids.
    outvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    tempvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
        # nothing needs to be changed in the next line for the neighbor problem.
        # NOTE - now StateRecord returns a Dict{Int64, NTuple{9, Int64}}
    stdict = StateRecord(D.all[el].nfids, el, D)                                      # returns the restricted state.
    for el2 in neighbors                                                              # adds keys for the neighbors to the temp dict. 
      outvals[D.all[el2].fid] = Dict{NTuple{10, Int64}, Float64}()
      tempvals[D.all[el2].fid] = Dict{NTuple{10, Int64},Float64}() 
      totest[D.all[el2].fid] = true                                                   # don't test convergence of neighbors temporarily.  FIXME 
      StateEnumerate(TupletoCNS(stdict[D.all[el2].fid]), outvals[D.all[el2].fid]) 
      StateEnumerate(TupletoCNS(stdict[D.all[el2].fid]), tempvals[D.all[el2].fid])
      locs[D.all[el2].fid] = el2
    end 
    StateEnumerate(D.all[el].cns, outvals[D.all[el].fid])                             # TODO - starting values here.
    StateEnumerate(D.all[el].cns, tempvals[D.all[el].fid])                            # this does NOT need starting values.  
    totest[D.all[el].fid] = true                                                      # all facilities to do initially set to false.  
    locs[D.all[el].fid] = el                                                          # stores a fid,location value
    push!(neighbors, el) # now this contains the location in D.all from chunk.
  end
  println(keys(locs))
  # Updating process:
  converge = false
  while (!converge)&(its<itlim)                                                       # if true keep going.    
    converge = true                                                                   # reassign, to catch when it terminates.
    for k in keys(totest)                                                             # TODO - not updating the competitors.
      if totest[k]                                                                    # only run those for which true.
        if messages println("tempvals keys before: ", keys(tempvals)) end
        if messages println("outvals keys before: ", keys(outvals)) end
        if messages println("locs keys before: ", keys(locs)) end 
        ExactChoice(tempvals, outvals, locs, k, locs[k], p1, p2, D; messages = true)  
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




#=


"""
`ExactVal(D::DynState, V::allvisits, itlim::Int64, chunk::Array{Int64,1}; debug::Bool = true)`
Computes the exact solution for smaller markets.  1 - 5 firms at most.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);

entries to consideR: 1-15 are all.
Here number of neighbors and entry (dyn.all[x])
Duopoly:
1
4
5
9

0 6 (monopoly)
0 7 (monopoly)
0 8 (monopoly)

triopoly:
18
20
21
24
25
26

4 13
4 14
4 15

#NB: consider an "itlim" ceiling

# testing: 

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);
ch = [1] # first element
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
ExactVal(dyn, ch, p1, p2)

PatientZero(p1, p2)

NOTES on current problems:
- some firms  are getting added... that is, fids are getting added which I don't want added.  Where does that happen?  

"""




=#
