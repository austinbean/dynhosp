# Now solve the dynamic model itself.


# Can put these into a market and then into a state.

# use the type neighbors !

using Distributions
using DataFrames
type hitcount # this type will record visits to a state-action pair
  conf::neighbors
  visits::Dict{Int64,Int64}
end

type history
  path::Dict{neighbors, hitcount}
  totalcount::Int64 # records total number of iterations
end


# Initialize empty:  hithist = history(Dict{neighbors, hitcount}(), 0)

type hstate
  #NB: track these for the level !
  conf::neighbors
  wvalues::Array{Float64,1}
  actprobs::Array{Float64,1}
end

type shortrec<:ProjectModule.Fac
  # needs to take location.  must measure distances.  need to update all of these every time, probably.
  fid::Int64
  lat::Float64
  long::Float64
  level::Int64
  truelevel::Int64
  ns::neighbors
  choices::Array{Int64, 2}
  chprobs::WeightVec
end


type simh<:ProjectModule.Fac
  fid::Int64
  level::Int64
  actual::Int64
  visited::Dict{neighbors, hstate} #possible given "isequal" and "hash" extended for "neighbors"
  ns::Array{shortrec, 1}
  exit::Bool
end


type DynState # this should hold a collection of ALL of the hospitals, for the approximation.
  all::Array{simh, 1}
end



"""
`DynStateCreate(Tex::EntireState)`
Create the dynamic records from the existing state, don't bother doing it from scratch.
And these don't have to be organized by zip.
Use the EntireState from the equilibrium simulation, not the first counterfactual.
Make using
TexasEq = MakeNew(ProjectModule.fips, ProjectModule.data05)
"""
function DynStateCreate( Tex::EntireState )
  outp = DynState(Array{simh,1}())
  for k1 in keys(Tex.mkts)
    println(k1)
    for hk in keys(Tex.mkts[k1].collection)
      newsimh = simh(Tex.mkts[k1].collection[hk].fid,
                     Tex.mkts[k1].collection[hk].level,
                     Tex.mkts[k1].collection[hk].level,
                     Dict{neighbors,hstate}(),
                     Array{shortrec,1}(),
                     false)
      for k2 in keys(Tex.mkts)
        for hk2 in keys(Tex.mkts[k2].collection)
          if hk2 != hk
            d1 = distance(Tex.mkts[k1].collection[hk].lat, Tex.mkts[k1].collection[hk].long , Tex.mkts[k2].collection[hk2].lat , Tex.mkts[k2].collection[hk2].long)
            if d1 < 25 #check distance
              push!(newsimh.ns, shortrec(Tex.mkts[k2].collection[hk2].fid,
                                         Tex.mkts[k2].collection[hk2].lat,
                                         Tex.mkts[k2].collection[hk2].long,
                                         Tex.mkts[k2].collection[hk2].level,
                                         Tex.mkts[k2].collection[hk2].level,
                                         neighbors(0,0,0,0,0,0,0,0,0),
                                         ChoicesAvailable(Tex.mkts[k2].collection[hk2]),
                                         Tex.mkts[k2].collection[hk2].chprobability))
            end
          end
        end
      end
      for el in newsimh.ns
        for el2 in newsimh.ns
          if el.fid != el2.fid
            if (distance(el.lat, el.long, el2.lat, el2.long)>0)&(distance(el.lat, el.long, el2.lat, el2.long)<5)
              if el2.level == 1
                el.ns.level105+=1
              elseif el2.level == 2
                  el.ns.level205+=1
              else # equals 3
                  el.ns.level305+=1
              end
            elseif (distance(el.lat, el.long, el2.lat, el2.long)>=5)&(distance(el.lat, el.long, el2.lat, el2.long)<15)
              if el2.level == 1
                  el.ns.level1515+=1
              elseif el2.level == 2
                  el.ns.level2515+=1
              else # equals 3
                  el.ns.level3515+=1
              end
            elseif (distance(el.lat, el.long, el2.lat, el2.long)>=15)&(distance(el.lat, el.long, el2.lat, el2.long)<=25)
              if el2.level == 1
                  el.ns.level11525+=1
              elseif el2.level == 2
                  el.ns.level21525+=1
              else # equals 3
                  el.ns.level31525+=1
              end
            end
          end
        end
      end
      push!(outp.all, newsimh)
    end
  end
  return outp
end



# Function to check the facilities with no neighbors:
"""
`NCheck(d::DynState, e::EntireState)`
This function just prints out the records of hospitals which have one or two neighbors to check them against
the same record for the same hospital in the EntireState record type.    Solely for debugging purposes.
"""
function NCheck(d::DynState, e::EntireState)
  for el in d.all
    if size(el.ns,1)<2
      println("*************************")
      println(el.fid)
      println(el.ns)
      println(e.mkts[e.fipsdirectory[el.fid]].collection[el.fid].hood)
    end
  end
end


"""
`GetProb`

"""


"""
`Simulate(d::DynState)`
This should take a dynamic state, as generated above in `DynStateCreate` and do a simulation.
What does that mean?
- Start in some state for each firm
- Choose actions by the other firms in the simh.ns
- Compute the demand and the return for the simulating hospital
- Write out the return to this, updating the estimate of the return at that state.
"""
function Simulate(d::DynState)



end




# The "Payoff" function in counter 1 gets the profit

"""
`ComputeR(hosp::simh, ppats::Dict{Int64, ProjectModule.patientcount}, mpats::Dict{Int64, ProjectModule.patientcount}, Tex::EntireState, wtp::Dict{Int64,Float64} )`
Computes the return (current profit + expected continuation) for each hospital in the state.
"""
function ComputeR(hosp::simh,
                  ppats::Dict{Int64, ProjectModule.patientcount},
                  mpats::Dict{Int64, ProjectModule.patientcount},
                  Tex::EntireState,
                  wtp::Dict{Int64,Float64} ;
                  disc::Float64 = 0.95)
  #NB: all of this can be updated in place, I think.
  for

    Payoff() + disc*(dot(hosp.wvalues, hosp.actprobs) + eulergamma - dot(log(hosp.actprobs),hosp.actprobs))

  end
end

"""
`PolicyUpdate(hosp::simh, neww::Array{Float64,1})`
Takes an array of the new W values and maps back to the simh action probabilities.
"""
function PolicyUpdate(hosp::simh, neww::Array{Float64,1})
  hosp.actprobs = exp(neww - maximum(neww))/sum(exp(neww-maximum(neww)))
end


"""
`UpdateW(hosp::simh, ret::Float64, stcounter::Int64)`
"""
function UpdateW(hosp::simh, ret::Float64, stcounter::Int64)


end



function CheckConvergence()


end
