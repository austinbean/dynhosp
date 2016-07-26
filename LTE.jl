# Laplace-type Estimator
# include("/Users/austinbean/Desktop/dynhosp/LTE.jl")


using DataFrames
using Distributions
using StatsBase

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


        function MetropolisHastings(initialpr::Vector,
                                    max_iterations::Int64,
                                    tolerance::Float64;
                                    param_dim = length(initialpr),
                                    pro_μ = zeros(param_dim),
                                    pro_σ = 10*eye(param_dim),
                                    proposal = Distributions.MvNormal(pro_μ, pro_σ),
                                    prior_μ = 500*ones(param_dim),
                                    prior_σ = 100*eye(param_dim),
                                    prior = Distributions.MvNormal(prior_μ, prior_σ))
          # Basics
          converged = false
          curr_it = 1

          # Storing the values:
          path = zeros(max_iterations, param_dim)

          # initial guess:
          curr_x = initialpr

          # Probability of initial guess according to prior
          curr_prior = pdf(prior, curr_x)

          # Value of objective function at initial guess
          curr_vals = objfun(curr_x)

          while curr_it < max_iterations && !converged # convergence flag not used yet
            # Proposed next value -
            # TODO: this is wrong: this must be *conditional*
            next_x = rand(proposal)

            # Probability of proposal at prior
            next_prior = pdf(prior, next_x)

            # Value of objective at Proposal
            # TODO: I think the problem is in the dimension of the output of the proposal distribution.
            next_vals = objfun(next_x)

            # Difference in value of objectives
            val_diff = exp(next_vals - curr_vals)

            rho = minimum( [val_diff*(next_prior/curr_prior)*(1) , 1.0]) # last term in parens 1 since dist chosen symmetric

            accept = sample([true false], WeightVec([rho, 1-rho]))

            if accept
              # Accepted Proposal
              curr_x = next_x
              curr_prior = next_prior
              curr_vals = next_vals
            end
            path[curr_it, :] = curr_x
            curr_it += 1
            # convergence flag not yet used.
          end

          return path

        end # of MetropolisHastings()

    sim_vals = MetropolisHastings(1000*ones(params), 100_000, 1e-12)

    # Results to return:
    println(mean(sim_vals,1))

end # of LTE function


LTE()
