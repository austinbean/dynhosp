



# Inspired by Optim.jl Simulated Annealing but with gradient information
# Several functions stolen from Optim.

macro sgtrace() #stochastic gradient trace
   quote
       if tracing
           dt = Dict()
           if o.extended_trace
               dt["x"] = copy(x)
               dt["g(x)"] = copy(g)
           end
           g_norm = vecnorm(g, Inf)
           update!(tr,
                   iteration,
                   f_x,
                   g_norm,
                   dt,
                   o.store_trace,
                   o.show_trace,
                   o.show_every,
                   o.callback)
       end
   end
end

#=
# These are specific functions
# Use this in the DifferentiableFunction construction.

function objfun(x::Vector; inp1::Array{Float64,2}=eq_opt, inp2::Array{Float64,2}=neq_opt)
  sum((min(inp1*x - inp2*x, 0)).^2)
end
=#
function objfun∇(x::Array{Float64,2}; inp1::Array{Float64,2}=eq_opt, inp2::Array{Float64,2}=neq_opt)
   grad = zeros(size(x))
   rows, columns = size(inp1)
     for k = 1:rows
       for j = 1:columns
         if ((inp1[k,j] - inp2[k,j])*x[j]>0)
           grad[j] += 2*(inp1[k,j] - inp2[k,j])*()
         end
       end
     end
   return grad
end


immutable StochasticGradient <: Optimizer
   direction::Function
end


# maybe this should go in "direction"
function neighbor!(x::Array, k::Int; grad::Function=objfun∇; Aparam::Int= 20, bparam::Float64=0.5, alpha::Float64=0.6, Bparam::Int= 20)
   x_proposal = zeros(size(x))
   # Implements a stochastic gradient procedure discussed in Spall (2003) - Stochastic Search and Optimization.
   a = grad*Aparam/(k+1+Aparam)^alpha # the gradient accessed as difffunction.g!
   b = (bparam/(sqrt(((k+1)^alpha)*log((k+1)^(1-alpha)+Bparam))))*randn(size(x))
   for i in 1:length(x)
        @inbounds x_proposal[i] = x[i] - a[i] + b[i]
   end
   return x_proposal
end

StochasticGradient(; direction::Function = neighbor!) =
StochasticGradient(neighbor!)

function optimize{T}(d::DifferentiableFunction,
                    initial_x::Vector{T},
                    mo::StochasticGradient,
                    o::OptimizationOptions)
   # initialize
   iterations = 0

   # Print header if show_trace is set
   print_header(o)

   # Maintain current state in x and previous state in x_previous
   x, x_previous = copy(initial_x), copy(initial_x)

   # Maintain current intermediate state in y and previous intermediate state in y_previous
   y, y_previous = copy(initial_x), copy(initial_x)

   # Count the total number of iterations
   iteration = 0

   # Track calls to function and gradient
   f_calls, g_calls = 0, 0

   # Count number of parameters
   n = length(x)

   # Maintain current gradient in g
   g = Array(T, n)

   # Store f(x) in f_x
   f_x_previous, f_x = NaN, d.fg!(x, g)
   f_calls, g_calls = f_calls + 1, g_calls + 1

   # Trace the history of states visited
   tr = OptimizationTrace(mo)
   tracing = o.store_trace || o.show_trace || o.extended_trace || o.callback != nothing
   @sgtrace

   # Assess types of convergence
   x_converged, f_converged, g_converged = false, false, false

   # Iterate until convergence
   converged = false

   # Call the trace for the first iteration
   @sgtrace()

 while !converged && iterations < max_iter



     iterations += 1
     # Gradient + random element:
     mo.direction( , iterations, d.g!(x,  ))

     @simd for i in 1:n
                 @inbounds x[i] = x[i] + dx[i]
     end

     # Assess Convergence:
     x_converged,
     f_converged,
     g_converged,
     converged = assess_convergence(x,
                                    x_previous,
                                    f_x,
                                    f_x_previous,
                                    g,
                                    o.x_tol,
                                    o.f_tol,
                                    o.g_tol)

     @sgtrace()

 end


return MultivariateOptimizationResults("Stochastic Gradient",
                                          initial_x,
                                          best_x,
                                          Float64(best_f_x),
                                          iteration,
                                          iteration == o.iterations,
                                          false,
                                          NaN,
                                          false,
                                          NaN,
                                          false,
                                          NaN,
                                          tr,
                                          f_calls,
                                           0)

end



###
