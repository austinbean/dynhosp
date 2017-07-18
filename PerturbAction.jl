# PerturbAction
#using Distributions
type VarianceException <: Exception end

function perturb(probs::Array, eps::Float64, control::Bool)
    # TODO - this is always four!
  actions = maximum(size(probs)) # always a column vector
  if control == true
    srand(123) #enables debugging by using seed.
  end
  if (eps > 1) | (eps < 0)
    return VarianceException
  end

  d = Normal(0, eps)
  # why transpose?
  pturb = rand(d, actions)'
  # TODO - this is a dumb thing to do. 
  return (max((probs + pturb),0))/(sum(max(probs + pturb, 0)))
end
