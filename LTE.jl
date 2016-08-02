# Laplace-type Estimator
# include("/Users/austinbean/Desktop/dynhosp/LTE.jl")


using DataFrames
using Distributions
using StatsBase
using Gadfly

function LTE()

    # Imports the simulation results.
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

    # Set other values to 0
    input1 = 0; interim_private_eq = 0; interim_private_neq = 0; interim_medicaid_eq = 0; interim_medicaid_neq = 0;
    eq_medicaid_lev1 = 0; eq_medicaid_lev2 = 0; eq_medicaid_lev3 = 0; neq_medicaid_lev1 = 0; neq_medicaid_lev2 = 0; neq_medicaid_lev3 = 0;


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
    plot(x=ones(20), Geom.point, Guide.xlabel("Die") )

    neq_opt = convert(Array{Float64, 2}, fout12[:, 4:end]);
    opt = eq_opt - neq_opt;
    hsims = 500 #size(fout1)[1] # number of simulations
    ncols = size(eq_opt,2) # number of columns with nonzeros (pairs!)
    const params = convert(Int, ncols) # don't think this conversion is strictly necessary


        function objfun(x::Vector; scale_fact = 1/10000, inp1::Array{Float64,2}=scale_fact*eq_opt, inp2::Array{Float64,2}=scale_fact*neq_opt, diffmat::Array{Float64,2}=inp1-inp2)
          sum(min(diffmat*x, 0).^2)
        end


        function MetropolisHastings(initialpr::Vector,
                                    max_iterations::Int64,
                                    tolerance::Float64;
                                    param_dim = length(initialpr),
                                    pro_μ = zeros(param_dim),
                                    pro_σ_scale::Float64 = 1.0,
                                    pro_σ = pro_σ_scale*eye(param_dim),
                                    proposal = Distributions.MvNormal(pro_μ, pro_σ),
                                    prior_μ_scale::Float64 = 1.0,
                                    prior_μ = prior_μ_scale*ones(param_dim),
                                    prior_σ_scale::Float64 = 1.0,
                                    prior_σ = prior_σ_scale*eye(param_dim),
                                    prior = Distributions.MvNormal(prior_μ, prior_σ))
# TODO: What's going on with the scaling of the parameters and the proposal distributions?  How do I
# rescale at the end?  And what is a proper distribution for the prior and the proposals?



#TODO: in the main loop below, am I getting the difference between the prior and proposal Distributions
# mixed up?  And should I include both in the log Hastings ratio? 

          # Basics
          converged = false
          curr_it = 1
          overflowcount = 0
          accepted = 0

          # Storing the values:
          path = zeros(max_iterations, param_dim) # 10 allocations / 50 MB

          # initial guess:
          curr_x = initialpr

          # Probability of initial guess according to prior
          curr_prior = pdf(Distributions.MvNormal(curr_x, prior_σ), curr_x) # 10 allocations.

          # Value of objective function at initial guess
          curr_vals = objfun(curr_x) # 14 allocations / 9 kb

          while curr_it < max_iterations && !converged # convergence flag not used yet
            # Proposed next value -
            next_x = curr_x + rand(proposal) #10 allocations / 970 bytes

            # Probability of proposal at prior
            next_prior = pdf(Distributions.MvNormal(curr_x, prior_σ), next_x) # 10 allocations / 1 kb

            # Value of objective at Proposal
            next_vals = objfun(next_x) #14 allocations / 9 kb

            # Difference in value of objectives
            #TODO: Consider the LOG of the hastings ratio here.  See Hdbk MCMC p.23
            val_diff = next_vals - curr_vals

            if val_diff == Inf  # keep track of times when this is too large.
              overflowcount += 1
            end

            # 8 allocations
            rho = minimum( [val_diff+log(next_prior)-log(curr_prior), 1.0])
            # This is accepting way too many.
        #    accept = sample([true false], WeightVec([rho, 1-rho])) #8 allocations

            if rand() < rho
              # Accepted Proposal
              curr_x = next_x
              curr_prior = next_prior
              curr_vals = next_vals
              accepted += 1
            end
             # Record the current parameter values whether they changed or not.
            for el = 1:param_dim
              @inbounds path[curr_it, el] = curr_x[el] # 4 allocations.
            end
            curr_it += 1
            # convergence flag not yet used.  What is it going to do?  Anything?
          end

          return path, overflowcount, accepted

        end # of MetropolisHastings()
    const nsims = 100000
    sim_vals, counter, accept = MetropolisHastings(ones(params), nsims, 1e-12)

    # Results to return:
    println(mean(sim_vals,1))
#    plot(x=rand(20), Geom.point, Guide.xlabel("This is a plot  "))
    println("Fraction Accepted ", accept/nsims)
    println("Count of Numerical Overflow ", counter)
  p1 = plot(x=sim_vals[:,1], Geom.histogram)
  p2 = plot(x=sim_vals[:,2], Geom.histogram)

  #=
  # for printing in the future:
  for el in 1:size(Optim.minimizer(result3), 1)
    print(varcolnames[el+3], "  ", Optim.minimizer(result3)[el], " param: ", paramsymbs[el], " symbol number: ", el, "\n")
  end
  =#
  return sim_vals

end # of LTE function


sims = LTE()

# function test()
#   x1 = rand(22, 130)
#   plot(x=mean(x1,1), Geom.histogram)
# end
# test()
