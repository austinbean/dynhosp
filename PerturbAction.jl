# PerturbAction
using Distributions

function perturb(probs::Array, eps::Float64, control::bool)
  actions = size(probs)[1]
  if control == true
    srand(123)
  elseif control == false 
    # nothing
  end

  d = Normal(0, eps)
  pturb = rand(d, actions)
  probs = probs + pturb

end
