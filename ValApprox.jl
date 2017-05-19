
"""
`ValApprox(D::DynState)`
This computes the dynamic simulation across all of the facilities in all of the markets.

TODO - where is the shock added and what is the variance of the shock?  It should be the normalizing constant.

To start:
dyn = CounterObjects(50);
V = allvisits(Dict{Int64, vrecord}());
ValApprox(dyn, V, 100_000 ; chunk = [2]) # just doing one hospital.
22.001426 seconds (281.50 M allocations: 5.099 GiB, 6.99% gc time)
"""
function ValApprox(D::DynState, V::allvisits, itlim::Int64; chunk::Array{Int64,1} = collect(1:size(D.all,1)), debug::Bool = false)
  iterations::Int64 = 0
  converged::Bool = false
  a::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  b::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  #steadylevs = AllAgg(D, chunk)
  for el in chunk                                                          # creates a dictionary of visited records.
    V.all[D.all[el].fid] = vrecord( Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, 1}(), 1)
  end
  while (iterations<itlim)&&(!converged)
    for el in D.all[chunk]
      if !el.converged                                                     # only keep simulating with the ones which haven't converged
        DSimNew(el.mk, el.fid, a, b)
        GetProb(el)                                                        # this chooses the action by the other firms
        # TODO: add a market size check here to do something more like oblivious in large markets.
        if !haskey(el.visited, KeyCreate(el.cns, el.level))
          println("adding new counter")
          el.visited[KeyCreate(el.cns, el.level)]=nlrec(MD(ChoicesAvailable(el), StartingVals(el, a, b)), vcat(ChoicesAvailable(el),transpose(PolicyUpdate(StartingVals(el, a, b)))), Dict(k => 0 for k in ChoicesAvailable(el)) )
        end
        Action = ChooseAction(el)                                              # Takes an action and returns it.
        ComputeR(el, a, b, Action, iterations; debug = debug)                  # Computes the return to the action
        level::Int64 = LevelFunction(el, Action)                                       # Level may change with action, but for next period.
        if iterations <= 1_000_000
          push!(V.all[el.fid].visited, RTuple(el, Action))                     # Record the first million state-action pairs in a vector
        elseif iterations >1_000_000
          V.all[el.fid].visited[iterations%1_000_000] = RTuple(el, Action)     # Once this is a million entries long, start overwriting to keep track of only 1_000_000
        end
        el.previous = el.level                                              # Reassign current level to previous.
        el.level = level                                                    # Reassign current level, if it has changed or not.
        ExCheck(el)                                                         # Checks for exit
        FixNN(el)                                                           # Fixes the firms neighbors.
        iterations += 1                                                     # Update iteration count - TODO: delete after debugging.
        V.all[el.fid].totalcnt += 1                                         # Update the iteration count within the visit records.
        PatientZero(a,b) # resets both patientcounts to zero.
      end
      # TODO - uncomment convergence test when that is debugged.
      if iterations%1_000_000 == 0                                        # Check for convergence every million iterations
        CheckConvergence(el) # FIXME - this is the wrong signature for this function.
      end
    end
  end
  converged = Halt(D, chunk)                                                # Check to see if all firms in "chunk" have converged, then halt if they have.
end
