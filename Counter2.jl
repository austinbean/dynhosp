# Now solve the dynamic model itself.


# Can put these into a market and then into a state.


type simh<:Fac
  level::Int64
  wvalues::Array{Float64,1}
  actprobs::Array{Float64,1}
  exit::Bool
end


# The "Payoff" function in counter 1 gets the profit

function ActionValue(hos::simh, pats::Dict{Int64, patientcount})



  return
end



function ActionProbs(hos::simh)


  return 
end


function CheckConvergence()


end
