"""
`ExactConvergence(current::Dict{ Int64, Dict{NTuple{10, Int64}, Dict{Int64, Float64} } }, stable::Dict{ Int64, Dict{NTuple{10, Int64}, Dict{Int64, Float64} } }; toler::Float64 =0.001, debug::Bool = true  )`
This will check convergence.  Does this by measuring the maximum difference at every state/action pair 
for each firm.  Returns a boolean recording convergence, but also returns a list of fids of unconverged facilities.
Operates on two dictionaries: one the permanent ("stable") and the other the temporary ("current")

Testing: 
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);;
test1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }();
test2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }();
test1[dyn.all[6].fid] = Dict{NTuple{10, Int64},  Float64 }();
test2[dyn.all[6].fid] = Dict{NTuple{10, Int64},  Float64 }();
StateEnumerate(dyn.all[6].cns, test1[dyn.all[6].fid])
StateEnumerate(dyn.all[6].cns, test2[dyn.all[6].fid])

test1[dyn.all[6].fid][(0,0,0,0,0,0,0,0,0,1)] = 20 #assign a value.
totest = Dict{Int64,Bool}()
totest[dyn.all[6].fid] = false 
ExactConvergence(test1, test2, totest; messages = false) == false # returns true.

ExactConvergence(test1, test2, totest; messages = true)
false # this is returning "converged" FALSE and the list of the unconverged facilities (in this case only one.)
"""
function ExactConvergence(current::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                          stable::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } },
                          totest::Dict{Int64,Bool}, # when the bool is "true", this has converged.
                          its::Int64; 
                          start_ch::Int64 = 10, # start checking convergence when iterations exceed this threshold.
                          messages::Bool = true,
                          toler::Float64 =1e-5)
  converge::Bool = false 
  diffs::Dict{Int64,Float64} = Dict{Int64,Float64}()                     # check only the guys still being done.
  # FIXME - would this problem be solved by initializing this dict with positive values?  
  for k1 in keys(stable)
    diffs[k1] = 1.0 # this is above toler
  end 
  if its > start_ch                                                      # provides a floor before which convergence is not checked and not updated. 
    for fid in keys(current)                                             # checks a subset ONLY, given by those in "current" whose locations are in chunk.  
      if !totest[fid]                                                    # keys in totest for which false (i.e., not converged)
        maxdiff::Float64 = -1.0   
        for state in keys(current[fid])                                  # states available to the firm.
          if haskey(stable[fid], state)
            if (current[fid][state]>0.0)&(stable[fid][state]>0.0)            # test states at which value is > 0, since some state values are not computed
              if current[fid][state] == stable[fid][state]
                # FIXME - too many of these are identical.  That can't be right.  
                # println("identical at: ", fid, " ", state)
              end 
              if abs(current[fid][state] - stable[fid][state]) > maxdiff # we want MAX difference. 
                maxdiff = abs(current[fid][state] - stable[fid][state])
              else 
                # do nothing - do not reassign diff.  
              end
            else 
              # do nothing, i.e., do not test when one or the other is 0.0
            end  
          else  
            if messages println("a state wasn't found ") end           # FIXME - eventually delete this branch.
          end 
        end
        #if messages println("for fid: ", fid) end
        if messages println("max diff: ", maxdiff) end
        diffs[fid] = abs(maxdiff)                                           # keep track of the max diff.  
      end 
    end 
    # FIXME - this can't be working right.  The next line is never printing.  
    if messages 
      println("current differences ")
      println(diffs)
    end 
    # Check for convergence - 
    # look at the maximum difference across states for each firm.
    # FIXME - now thre problem is: everything gets initialized to zero, convergence achieved right away, never tested again. 
    for k1 in keys(diffs)                                              # by iterating over keys in diffs, we only check facs which have not converged yet.  
      converge = converge&(diffs[k1]<toler)
      if diffs[k1] > toler                                             # totest exists in the scope of ExactVal, recording which firms are around and whether they have converged.  
        totest[k1] = false 
      else 
        totest[k1] = true 
      end 
    end 
  end # end of "if its" condition.
end 



#=

test1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }();
test2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }();
test1[4450450] = Dict{NTuple{10, Int64},  Float64 }();
test2[4450450] = Dict{NTuple{10, Int64},  Float64 }();

for i = 0:9
  test1[4450450][(i,i,i,i,i,i,i,i,i,i)] = rand()
  test2[4450450][(i,i,i,i,i,i,i,i,i,i)] = 0.0
end 


function DCopTest(d1, d2)
  for k1 in keys(d1)
    println("TEST1:")
    for k2 in keys(d1[k1])
      println(k2, " => ", d1[k1][k2])
      d2[k1][k2] = d1[k1][k2]
      d1[k1][k2] = 0
    end
  end  
end 

DCopTest(test1, test2)
test1[4450450]
test2[4450450]


=#




