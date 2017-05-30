
"""
`ValApprox(D::DynState)`
This computes the dynamic simulation across all of the facilities in all of the markets.
- Results are stored in the Dyn record.  That is, values.
- These are accessed as dyn.ll[#].visited.[(neighbors, level)].aw
- the probabilities associated to the values are in dyn.all[#].visited[(neighbors, leve)].psi 
- there is a count of each state in dyn.all[#].visited[(neighbors,level)].visited.
- V is of type allvisits, which records which states was visited in the last million iterations.  

TODO - where is the shock added and what is the variance of the shock?  It should be the normalizing constant.

TODO - probabilities are not getting recorded properly: 

dyn.all[1].visited[(0,0,0,0,0,0,1,0,0,1)].aw
Dict{Int64,Float64} with 4 entries:
  10 => 0.0228582
  2  => 0.0228582
  11 => 2.28582e-5
  1  => 0.0228582

BUT:

  dyn.all[1].visited[(0,0,0,0,0,0,1,0,0,1)].psi
2×4 Array{Float64,2}:
 10.0       2.0       1.0       11.0
  0.251419  0.251419  0.251419   0.245743

TODO - counter is not getting incremented.  

dtest = Dict(1 => 0)

function testd(n::Int64, d1::Dict)
  for i = 1:n 
    for k1 in keys(d1)
      for j = 1:10
        d1[k1] += 1
      end 
    end 
  end 
end 

testd(100, dtest)


To start:
dyn = CounterObjects(50);
V = allvisits(Dict{Int64, vrecord}());
ValApprox(dyn, V, 1000 ; chunk = [2]) # just doing one hospital.
22.001426 seconds (281.50 M allocations: 5.099 GiB, 6.99% gc time)
"""
function ValApprox(D::DynState, V::allvisits, itlim::Int64; chunk::Array{Int64,1} = collect(1:size(D.all,1)), debug::Bool = false)
  iterations::Int64 = 0
  converged::Bool = false
  a::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  b::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  for el in chunk                                                          # creates a dictionary of visited records.
    V.all[D.all[el].fid] = vrecord( Array{NTuple{11,Int64}, 1}(), 1)
  end
  while (iterations<itlim)&&(!converged)
    for el in D.all[chunk]
      if !el.converged                                                     # only keep simulating with the ones which haven't converged
        DSimNew(el.mk, el.fid, a, b)
        GetProb(el)                                                        # this chooses the action by the other firms
        if !haskey(el.visited, KeyCreate(el.cns, el.level))
          println("adding new counter")
          el.visited[KeyCreate(el.cns, el.level)]=nlrec(MD(ChoicesAvailable(el), StartingVals(el, a, b)), PolicyUp2(ChoicesAvailable(el),PolicyUpdate(StartingVals(el, a, b)))), Dict(k => 1 for k in ChoicesAvailable(el)) )
        end
        Action = ChooseAction(el)                                              # Takes an action and returns it.
        if Action!=10
          println("action: ", Action )
        end 
        # TODO: print out the keys in visited here and afterwards
        # DELETEBELOW 
        # for k1 in keys(el.visited)
        #   for k2 in keys(el.visited[k1].counter)
        #     println("BEFORE: ", k1, " ", k2, " ", el.visited[k1].counter[k2])
        #   end 
        # end 
        #DELETEABOVE
        ComputeR(el, a, b, Action, iterations; debug = debug)                  # Computes the return to the action
        #DELETEBELOW
        for k1 in keys(el.visited)
          for k2 in keys(el.visited[k1].counter)
            println("INSIDE: ", k1, " ", k2, " ", el.visited[k1].counter[k2])
          end 
        end 
        #DELETEABOVE 
        level::Int64 = LevelFunction(el, Action)                               # Level may change with action, but for next period.
        if iterations <= 1_000_000
          push!(V.all[el.fid].visited, RTuple(el, Action))                     # Record the first million state-action pairs in a vector
        elseif iterations >1_000_000
          V.all[el.fid].visited[iterations%1_000_000] = RTuple(el, Action)     # Once this is a million entries long, start overwriting to keep track of only 1_000_000
        end
        if level != el.level                                                   # levels don't agree - i.e., "level" here is the next level which has been drawn.
          println("level changed!")  
          UpdateDUtil(el)                                                      # this should update the utility for hospitals which changed level.  
        end 
        # FIXME - where is the probability being updated?  The function ProbUpdate - fix that. 
        # FIXME - where is the counter getting updated?  This was previously working.
        el.previous = el.level                                                 # Reassign current level to previous.
        el.level = level                                                       # Reassign current level, if it has changed or not.
        ExCheck(el)                                                            # Checks for exit
        FixNN(el)                                                              # Fixes the firms neighbors.
        iterations += 1                                                        # Update iteration count 
        V.all[el.fid].totalcnt += 1                                            # Update the iteration count within the visit records.
        PatientZero(a,b)                                                       # resets both patientcounts to zero.
      end
      if iterations%1_00 == 0                                                  # Check for convergence every million iterations
        CheckConvergence(el, V.all[el.fid].visited; debug = true) 
      end
    end
    #DELETEBELOW
    for el in D.all[chunk]
      for k1 in keys(el.visited)
        for k2 in keys(el.visited[k1].counter)
          println("AFTER MAIN: ", k1, " ", k2, " ", el.visited[k1].counter[k2])
        end 
      end 
    end 
    #DELETEABOVE
  end 
  # DELETE BELOW 
  for el in D.all[chunk]
    for k1 in keys(el.visited)
      for k2 in keys(el.visited[k1].counter)
        println("EXIT: ", k1, " ", k2, " ", el.visited[k1].counter[k2])
      end 
    end 
  end 
  #DELETE ABOVE
  converged = Halt(D, chunk)                                                # Check to see if all firms in "chunk" have converged, then halt if they have.
end
