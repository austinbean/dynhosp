# Computes the value of the simulation path, taking the state history as
# given.


function DynamicValue(state_history::Array, beta::Float64)
  len, width = size(state_history)
  fac_section = state_history[:, width-4]
  agg_section = state_history[:, width-3:end]

  mtotal = maximum(sum(agg_section[:,1:3 ], 2)) #maximum total hospitals.
  mlev1 = maximum(agg_section[:,1]) # maximum level 1
  mlev2 = maximum(agg_section[:,2]) # maximum level 2
  mlev3 = maximum(agg_section[:,3]) # maximum level 3


end
