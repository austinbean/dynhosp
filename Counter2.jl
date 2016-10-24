# Now solve the dynamic model itself.


# Can put these into a market and then into a state.

# use the type neighbors !

type state
  conf::neighbors
  wvalues::Array{Float64,1}
  actprobs::Array{Float64,1}
end

#=
function Base.isequal(n1::neighbors, n2::neighbors)::Bool
  n1.level105 == n2.level105 && n1.level205 == n2.level205 && n1.level305 == n2.level305 && n1.level1515 == n2.level1515 && n1.level2515 == n2.level2515 && n1.level3515 == n2.level3515 && n1.level11525 == n2.level11525 && n1.level21525 == n2.level21525 && n1.level31525 == n2.level31525
end

function Base.hash(n1::neighbors)
  hash((n1.level105, n1.level205, n1.level305, n1.level1515, n1.level2515, n1.level3515, n1.level11525, n1.level21525, n1.level31525) )
end
=#

type simh<:Fac
  level::Int64
  actual::Int64
  visited::Dict{neighbors, state}
  exit::Bool
end

#NB: create the state with the values.  Need good starting values

DynTex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())

function DynStateCreate( ; fi = ProjectModule.fips)


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
