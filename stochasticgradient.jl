# StochasticGradient
# Inspired by Optim.jl Simulated Annealing but with gradient information

#TODO: Parameters for the multiplicative elements A, B, etc.

function objfun(x::Vector; inp1::Array{Float64,2}=eq_opt, inp2::Array{Float64,2}=neq_opt)
   sum((min(inp1*x - inp2*x, 0)).^2)
end

function objfun∇(x::Array{Float64,2}; inp1::Array{Float64,2}=eq_opt, inp2::Array{Float64,2}=neq_opt)
    grad = zeros(size(x))
    rows, columns = size(inp1)
      for k = 1:rows
        for j = 1:columns
          if ((inp1[k,j] - inp2[k,j])*x[j]>0)
            grad[j] += 2*(inp1[k,j] - inp2[k,j])
          end
        end
      end
    return grad
end

function neighbor(x::Array, k::Int; Aparam::Int= 20, bparam::Float64=0.5, alpha::Float64=0.6, Bparam::Int= 20)
    x_proposal = zeros(size(x))
    # Implements a stochastic gradient procedure discussed in Spall (2003) - Stochastic Search and Optimization.
    a = objfun∇(x)*Aparam/(k+1+Aparam)^alpha
    b = (bparam/(sqrt(((k+1)^alpha)*log((k+1)^(1-alpha)+Bparam))))*randn(size(x))
    for i in 1:length(x)
         x_proposal[i] = x[i] - a[i] + b[i] #@inbounds
    end
    return x_proposal
end

while !converged

    converged = assess_convergence()

end










###
