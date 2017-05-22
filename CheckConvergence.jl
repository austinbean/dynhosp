

"""
`CheckConvergence(h::simh, V::Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},1}; draws::Int64 = 100, demands::Int64 = 10, disc::Float64 = 0.95, debug::Bool = true)`
Check the convergence criterion in Collard-Wexler or Pakes McGuire.
This can be checked every million iterations.  When that happens,
reset the counters for every state in "visited".
"Visited" accessed by V.all[fid].visited
- Convergence should be assessed at *visited* states, something like unique(V), not all of them.
"""
function CheckConvergence(h::simh, V::Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},1}; draws::Int64 = 1000, disc::Float64 = 0.95, debug::Bool = true)
  outp::Array{Float64,1} = Array{Float64, 1}()
  pairs::Array{Tuple{Float64,Float64},1} = Array{Tuple{Float64, Float64},1}()
  totvisits::Int64 = 0
  a::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  b::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  states::Dict{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, Tuple{Float64,Float64}} = Dict{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, Tuple{Float64,Float64}}()
  itercount::Int64 = 0
  for k in unique(V)                                                                                              # only check this set of values visited in the last million iterations.
   k1::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, k2::Int64 = KeytoTuple(k)
   approxim::Float64 = 0.0
   for d = 1:draws
      origlevel::Int64 = h.level                                                                          # keep track of the level inside of the loop so that it can be reset.
      nextact::Int64 = convert(Int64, sample(h.visited[k1].psi[1,:], ProjectModule.WeightVec(h.visited[k1].psi[2,:])))  # Take an action.  NB: LevelFunction takes Int64 argument in second place.
      h.level = LevelFunction(h, nextact)                                                                 # this level must be updated so that the profit computation is correct.
      DSimNew(h.mk, h.fid, a, b)
      currpi::Float64 = SinglePay(h, a, b, nextact)                                              # Current period return, excluding continuation value.
      contval::Float64 = 0.0
      if haskey(h.visited, KeyCreate(h.cns, h.level))                                                     # check neighbors/level pair
        #FIXME - note here: ContError should not be present upon exit.
        contval = disc*WProb(h.visited[KeyCreate(h.cns, h.level)])
        if nextact != 11
          contval += disc*(ContError(h.visited[KeyCreate(h.cns, h.level)]))
        end
      else
        #FIXME - what is happening here now on "not available"  contval stays 0 - but that's going to make the error larger.
        # But these will still be 0 since that's what StartingVals is giving.
        println("Added entry")
        # FIXME - this is a terrible line.  The closure will be really slow.  
        # but this only needs to occur when there is a state not seen before which is generated in the convergence check. 
        h.visited[KeyCreate(h.cns, h.level)]=nlrec(MD(ChoicesAvailable(h), StartingVals(h, a, b)), vcat(ChoicesAvailable(h),transpose(PolicyUpdate(StartingVals(h, a, b)))), Dict(k => 0 for k in ChoicesAvailable(h)) )
        if nextact!=11
          contval += 0 #disc*() #FIXME - this is not done.  
        end
      end
      h.level = origlevel                                                                                 # reset the level to the original value.
      approxim += (currpi+contval)                                                                        # this needs to be weighted by the right count
      PatientZero(a,b) # resets both patientcounts to zero.
   end
   push!(outp, (approxim/draws - h.visited[k1].aw[k2])^2)                                                 # TODO - replace this with a sum when confidence is reached in the outcome..
   push!(pairs, (approxim/draws, h.visited[k1].aw[k2]))
   states[KeyCreate(h.cns, h.level)] = (approxim/draws, h.visited[k1].aw[k2])
  end
  if debug
    for k1 in keys(h.visited)
      for k2 in keys(h.visited[k1].counter)
        itercount += h.visited[k1].counter[k2]                                                                 # how many iterations were made in total?
        h.visited[k1].counter[k2] = 1                                                                          # Reset all counter keys to 1 to track which states visited in last million iterations.
      end
    end
  end
  return outp, itercount, states, pairs                                                                             # TODO: eventually divide former by latter. And drop states.
end


