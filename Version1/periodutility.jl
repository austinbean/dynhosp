# period utility functions

function util(demand::Array, current::Array, previous::Array, thetas::Array, entry::Array, α2 = 0.1, α3 = 0.15 )
  # Array records entry at either lev 2 or 3 as [i, j] with i,j ∈ {0,1}
  # α_k, k ∈ {2,3} gives the fraction getting sent to level k facilities
  # Current is an array [i, j, k] ∈ {0,1}^3 recording current level, i+j+k = 1
  # Previous is an array like current recording what the firm did last period
  # Demand is an array α_k

# 0/1 recording current levels
  clev1 = current[1]
  clev2 = current[2]
  clev3 = current[3]

# 0/1 recording previous level
  plev1 = previous[1] # this isn't going to be used
  plev2 = previous[2]
  plev3 = previous[3]

# Theta parameter values:
  Θ_1 = thetas[1]
  Θ_2 = thetas[2]
  Θ_3 = thetas[3]

# entry:
  Entry2 = entry[1]
  Entry3 = entry[2]

# I need to make some assumptions about how demand will actually operate, but
# for now suppose that it comes in three levels.

  demand = clev1*demand[1] + clev2*demand[2] + clev3*demand[3]

# Revenue:
  rev = demand*(Θ_1 + clev2*α2*Θ_2 + clev3*α3*Θ_3)

# Entry costs:
  ent_costs = (clev2)*(1 - plev2)*Entry2 + clev3*(1-plev3)*Entry3


  return rev - ent_costs

end
