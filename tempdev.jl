function ExactVal(D::DynState,
                  chunk::Array{Int64,1}, 
                  p1::patientcount,
                  p2::patientcount;
                  messages::Bool = false,
                  itlim::Int64 = 350,
                  debug::Bool = true,
                  beta::Float64 = 0.95,
                  conv::Float64 = 0.0001)
  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  tempvals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  totest::Dict{Int64,Bool} = Dict{Int64,Bool}()                                       # will record convergence 
  all_locs::Dict{Int64,Int64} = Dict{Int64, Int64}()                                  # will record the locations of competitors
  st_dict::Dict{Int64,NTuple{9,Int64}} = Dict{Int64,NTuple{9,Int64}}()                # will record the states of all firms from the point of view of el.
  its::Int64 = 0                                                                      # records iterations, but will be dropped after debugging.
  for el in chunk # goal of this loop is to: set up the dictionaries containing values with entries for the fids.  
    neighbors::Array{Int64,1} = FindComps(D, D.all[el])                               #  these are addresses of fids.
    push!(neighbors, el)                                                              # add the location of the firm in chunk
    outvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    tempvals[D.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    CompsDict(neighbors, D, all_locs)                                                 # now this is Dict{Fid, Location}
    StateRecord(all_locs, D, st_dict)                                                 # returns the restricted state.
    for el2 in neighbors                                                              # adds keys for the neighbors to the temp dict. 
      outvals[D.all[el2].fid] = Dict{NTuple{10, Int64}, Float64}()
      tempvals[D.all[el2].fid] = Dict{NTuple{10, Int64},Float64}() 
      totest[D.all[el2].fid] = false                                                   # don't test convergence of neighbors temporarily.  FIXME 
      StateEnumerate(TupletoCNS(st_dict[D.all[el2].fid]), outvals[D.all[el2].fid]) 
      StateEnumerate(TupletoCNS(st_dict[D.all[el2].fid]), tempvals[D.all[el2].fid])
      #locs[D.all[el2].fid] = el2
    end 
    StateEnumerate(D.all[el].cns, outvals[D.all[el].fid])                              # TODO - starting values here.
    StateEnumerate(D.all[el].cns, tempvals[D.all[el].fid])                             # this does NOT need starting values.  
    totest[D.all[el].fid] = false                                                       # all facilities to do initially set to false.  
    #locs[D.all[el].fid] = el                                                          # stores a fid,location value
  end
  # Updating process:
  # NOTE - totest contains only the firms in the "market" which we want. 
  # FIXME - something needs to be done to initialize this...  
  converge::Bool = true
  while (converge)&&(its<itlim)                                                        # if true keep going.  ]
    if (its>328)&(its<335) # an error repeatedly happens after 330 iterations.
        # basically print everything...
        println("iteration: ", its, " fid ") 
        for k1 in keys(outvals)
            for k2 in keys(outvals[k1])
                if isnan(outvals[k1][k2])||isnan(tempvals[k1][k2])
                    # somehow this is happening after... 335 iterations?  Twice in a row?  
                    # five times in a row.  Always starting with 3490795.  This is weird.  
                    println(its, " FID: ",k1, " STATE: ", k2)
                end 
            end 
        end
        DSimNew( D.all[].mk, 3490795, p1, p2)
        println("Demand Sim: ", p1, "  ", p2)
        
    end    
#    converge = true                                                                    # reassign, to catch when it terminates.
    for k in keys(totest)                                                              # TODO - not updating the competitors.
      if totest[k]                                                                     # only run those for which true. 
        ExactChoice(tempvals, outvals, all_locs, st_dict, k, all_locs[k], p1, p2, D; messages = false)  
      end 
    end
    # if isnan(sum(p1))||isnan(sum(p2))
    #     println("from ExactVal: ")
    #     println("p1, p2: ", p1, "  ", p2)
    # end 
     
    # Convergence Test:
    # FIXME - I'll bet the problem is that it checks ALL states, some of which are initialized to zero, so 
    # these show as converged.  But could that matter?  It's supposed to focus on the max difference.  
    converge = ExactConvergence(tempvals, outvals, totest; messages = false)    
    # Copy the values and clean up.
    DictCopy(outvals, tempvals)
    DictClean(tempvals)                                                               # sets up for rewriting.
    for ky1 in keys(totest)                                                           # this tests every facility every time, but that's ok. 
      converge = ConvTest(totest)                                                     # iterates over bools in totest returns product
    end   
    its += 1
    # println("iteration ", its)
    # println("converge? ", converge)
  end 
  # Return equilibrium values:
  return outvals
end




#=


"""
`ExactVal(D::DynState, V::allvisits, itlim::Int64, chunk::Array{Int64,1}; debug::Bool = true)`
Computes the exact solution for smaller markets.  1 - 5 firms at most.


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


ch2 = [11]
ExactVal(dyn, ch2, p1, p2)

NOTES on current problems:
- some firms  are getting added... that is, fids are getting added which I don't want added.  Where does that happen?  

"""




=#
