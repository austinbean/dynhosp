# StochasticGradient
#=
Benchmark - Simulated Annealing (100_000 iterations)

9.635826 seconds (2.68 M allocations: 945.798 MB, 1.00% gc time)

100_000 iterations in this file:
(Some of the cost will be the file loading, but most probably not. )
389.816312 seconds (5.73 G allocations: 412.663 GB, 15.41% gc time)

function objfun∇(x::Vector)
  return ones(size(x)) #rewrite this.
end
=#

using DataFrames
using Distributions

type BestVal
  x::Array{Float64}
  val::Float64
end


function importstuff()

  input1 = readtable("/Users/austinbean/Google Drive/Simulation Results/combinedresults.csv");

  # This keeps the privately insured patient counts:
  # The division between first and second half is around 172.
    interim_private_eq = convert(Matrix, hcat(input1[:,1:81], input1[:,160:171]));
    interim_private_neq = convert(Matrix, hcat(input1[:,1:3], input1[:,172:249], input1[:,328:339]));
  # This keeps the Medicaid patient counts:
    interim_medicaid_eq = hcat(input1[:,1:3], input1[:,82:159]);
    interim_medicaid_neq = hcat(input1[:,1:3], input1[:,250:327]);
  # This keeps the equilibrium sim medicaid patient counts by level - interim_medicaid summed across competitor numbers
    eq_medicaid_lev1 = hcat(convert(Matrix, interim_medicaid_eq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_eq[:,4:29]),2));
    eq_medicaid_lev2 = hcat(convert(Matrix, interim_medicaid_eq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_eq[:,30:55]),2));
    eq_medicaid_lev3 = hcat(convert(Matrix, interim_medicaid_eq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_eq[:,56:81]),2));
  # This keeps the non-equilibrium sim medicaid counts by level
    neq_medicaid_lev1 = hcat(convert(Matrix, interim_medicaid_neq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_neq[:,4:29]),2));
    neq_medicaid_lev2 = hcat(convert(Matrix, interim_medicaid_neq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_neq[:,30:55]),2));
    neq_medicaid_lev3 = hcat(convert(Matrix, interim_medicaid_neq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_neq[:,56:81]),2));


  # The form of the output is the same as the input, except at each level the 26 medicaid columns have been summed into one
  # This creates TWO matrices, one with all of the equilibrium results, the other with all of the non-equilibrium results.
  fout11 = hcat( interim_private_eq[:,1:29], eq_medicaid_lev1[:,4], interim_private_eq[:,30:55], eq_medicaid_lev2[:,4], interim_private_eq[:,56:81], eq_medicaid_lev3[:,4], interim_private_eq[:,82:end])
  fout12 = hcat(interim_private_neq[:,1:29], neq_medicaid_lev1[:,4], interim_private_neq[:,30:55], neq_medicaid_lev2[:,4], interim_private_neq[:,56:81], neq_medicaid_lev3[:,4], interim_private_neq[:,82:end])

  varcolnames = Array{Symbol}(:0)
  push!(varcolnames, :fipscode)
  push!(varcolnames, :fid)
  push!(varcolnames, :year)
  for elem in ["EQ", "NEQ"]
    for name in ["PI"]
      for j = 1:3
        for k = 0:25
          push!(varcolnames, parse("$name"*"$elem"*"Lev$j"*"Comp$k"))
        end
          push!(varcolnames, parse("$elem"*"Medicaid"*"Lev$j"))
      end
    end
    for x in [1 2 3]
      for y in [1 2 3 "EX"]
        if x != y
          push!(varcolnames, parse("$elem"*"Trans$x$y"))
        end
      end
    end
    push!(varcolnames, parse("$elem"*"Enter1"))
    push!(varcolnames, parse("$elem"*"Enter2"))
    push!(varcolnames, parse("$elem"*"Enter3"))
  end

  #names!(fout1, colnames)

  paramsymbs = Array{UTF8String}(0)
  for k = 1:3
    for i = 0:25
      push!(paramsymbs, "Θ"*"$k"*"C$i")
    end
    push!(paramsymbs, "μ"*"$k")
  end
  for y in [1 2 3]
    for z in [1 2 3 "EX"]
      if y != z
        push!(paramsymbs, "ψ"*"$y"*"$z")
      end
    end
  end

  push!(paramsymbs, "Γ1")
  push!(paramsymbs, "Γ2")
  push!(paramsymbs, "Γ3")


  #=

  # This just checks that varcolnames and parasymbs are correct

  for i in 1:size(paramsymbs,1)
    println(varcolnames[3+i],"  ", paramsymbs[i])
  end

  # Delete columns of zeros - this is just for the testing part.  Eventually hopefully all will be filled in.
  # This only deletes if both the column for the equilibrium AND non-equilibrium are zero.
  # this doesn't work because the size is changing dynamically.
  # This will check quickly if there are pairs of columns which are all zeros.

  =#

  deletdsym = Array{Int64}(0) # vector of ints, for symbol indices (no "push!" method for arrays of symbols?)
  deletdcols = Array{Int64,1}(0) # vector of ints, for column indices - can index columns with this.  Start with ALL ints and remove these.
  pcop = deepcopy(paramsymbs)
  for col in 4:size(fout11,2)
    el = varcolnames[col] # these two lines give column names
    nel = varcolnames[col+93]
    if (sum(fout11[:,col]) == 0) & (sum(fout12[:,col]) == 0)
   #    println("Empty Column: ", col, "  ",sum(fout11[:,col]), "  ", el, "   ",sum(fout12[:,col]), "   ", nel )
    #  fout1 = fout1[:, (1:size(fout1,2).!=col)&(1:size(fout1,2).!=col+93)]
      push!(deletdcols, col)
      push!(deletdsym, col)
    else
    #  println(sum(fout11[:,col]), "  ",sum(fout12[:,col]) )
    end
  end
  deleteat!(paramsymbs, [x-3 for x in deletdsym]) # remember that the matrices include a 3-vector of identifiers.
  deleteat!(varcolnames, deletdsym) # removes the names of the colums - useful after the optimization to print the vals
  # Deleting columns:
  nonzers = setdiff(collect(1:size(fout11,2)), deletdcols) # keep only the columns not in the set deletdcols
  fout11 = fout11[:,nonzers]
  fout12 = fout12[:,nonzers]
  #=
  # This prints a list of the remaining coefficients, plus their indices in paramsymbs.
  for el in 1:size(pcop, 1)
    index = findfirst(paramsymbs, pcop[el])
    if index != 0
      println(index, "   ", paramsymbs[index], "   ", pcop[el])
    else
      println("*****   ", pcop[el])
    end
  end
  =#
  # Drop identifiers:
  eq_opt = convert(Array{Float64, 2}, fout11[:,4:end]);
  neq_opt = convert(Array{Float64, 2}, fout12[:, 4:end]);
  opt = eq_opt - neq_opt;
  hsims = 500 #size(fout1)[1] # number of simulations
  ncols = size(eq_opt,2) # number of columns with nonzeros (pairs!)
  const params = convert(Int, ncols) # don't think this conversion is strictly necessary


    function objfun(x::Vector; inp1::Array{Float64,2}=eq_opt, inp2::Array{Float64,2}=neq_opt)
       sum((min(inp1*x - inp2*x, 0)).^2)
    end

    function objfun∇(x::Array; inp1::Array{Float64,2}=eq_opt, inp2::Array{Float64,2}=neq_opt)
      params = length(x)
      gradient = zeros(params)
      diff = inp1 - inp2
      for i = 1:params # indexes columns
        for j = 1:size(diff,1) #indexes rows/firms
          if diff[j,i] < 0 # could be replaced in the next line with min(inp1[j,i] - inp2[j,i], 0) at the end
            gradient[i] += (2*sum(min(diff[j,:]*x, 0))*(diff[j,i])) # sum turns a 1x1 array into a scalar
          end
        end
      end
      return gradient
    end

#=
# Easy functions for testing.
    function objfun(x::Vector)
      return sum(x.^2)
    end

    function objfun∇(x::Vector)
      return 2*x
    end
=#

    function optimizeit(trace::Bool, iters::Int64; def_sd = 10, paramsize = params)
      # Parameters for the search coefficients
      parameter_dim = paramsize #some constant.

      # Tolerance:
      tol = 1e-10

      # Parameters for the coefficients in the stochastic search
      # Aparam::Int= 20;
      # bparam::Float64=0.5;
      # alp::Float64=0.6;
      # Bparam::Int= 20
      sd = def_sd # If you want the standard deviation of the normal disturbance to be larger.

      # Maximum Iterations allowed
      max_iter = iters

      #Initialize the count of iterations
      iterations = 1 # start at 1 since we take the log of this below.

      # Initial Guess for parameters, Initial Value for "next"
      x = 100*ones(parameter_dim)
      x_next = zeros(parameter_dim)

      # Initial value for function:
      f_previous = objfun(x)
      f_current = 0

      # Initial converged = false
      f_converged = false
      x_converged = false

      # Best values initialized at initial values.
      best = BestVal(x, objfun(x))

      while !f_converged && iterations < max_iter
        # an alternate set of coefficients
        # (Aparam/(iterations+1+Aparam)^alp)
        # (bparam/(sqrt(((iterations+1)^alp)*log((iterations+1)^(1-alp)+Bparam))))

        # Search in the direction
        ai = 1/(10*log(1+iterations))
        bi = 2/(1+iterations)
        x_next = x .- (ai/bi).*objfun∇(x).*rand(size(x))
        f_current = objfun(x_next)

        if trace
          if iterations % (max_iter/10) == 0
            println("Current iteration ", iterations)
            println("Mean current ", mean(x))
            println("Current x ", x_next)
          end
        end

        # Is current value better than best value?
        if f_current < best.val
          best.x = x_next
          best.val = f_current
        end

        # Has function converged?
        if abs(f_current - f_previous) < tol
          f_converged = true
        end
        # Has x converged?  This should almost never happen due to randomization - does not stop evaluation.
        if norm(x- x_next,1) < 1e-5
          x_converged = true
        end

        # assign current function value to "previous" for comparison next round
        f_previous = f_current
        # Assign current values of x to "x" for use next round
        x = x_next
        iterations += 1
      end
      println("Current Parameters ", x)
      println("Current Function Value ", objfun(x))
      println("Best Parameter Values ", best.x)
      println("Best Function Value ", best.val)
      println("Current Iteration: ", iterations)
      println("Reached Max Iterations: ", iterations == max_iter)
      println("Converged Flag: ", f_converged)
      println("X Converged Flag: ", x_converged)
    end

    return optimizeit(true, 10000; paramsize = params)


end # of main function.

importstuff()




###
