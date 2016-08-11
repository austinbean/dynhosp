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


        function objfun(x::Vector; scale_fact = 1/10000, inp1::Array{Float64,2}=scale_fact*eq_opt, inp2::Array{Float64,2}=scale_fact*neq_opt, diffmat::Array{Float64,2}=inp1-inp2)
          sum(min(diffmat*x, 0).^2)
        end


        function MetropolisHastings(initialpr::Vector,
                                    max_iterations::Int64;
                                    param_dim = length(initialpr),
                                    pro_μ = zeros(param_dim),
                                    pro_σ_scale::Float64 = 100.0,
                                    pro_σ = pro_σ_scale*eye(param_dim),
                                    proposal = Distributions.MvNormal(pro_μ, pro_σ),
                                    prior_μ_scale::Float64 = 1000.0,
                                    prior_μ = prior_μ_scale*ones(param_dim),
                                    prior_σ_scale::Float64 = 1000.0,
                                    prior_σ = prior_σ_scale*eye(param_dim),
                                    prior = Distributions.MvNormal(prior_μ, prior_σ),
                                    debug::Bool = true)
          # Basics
          curr_it = 2
          overflowcount = 0
          underflowcount = 0
          priorzerocount = 0
          accepted = 0
          if debug
            trace = zeros(max_iterations*param_dim, 6+param_dim)
            allvals = zeros(max_iterations*param_dim, param_dim)
            param_accept = zeros(param_dim)
            counter = 1
          end

          # Storing the values:
          path = zeros(max_iterations, param_dim) # 10 allocations / 50 MB
          for j = 1:param_dim
            path[1,j] = initialpr[j] #record initial guess
            allvals[1,j] = initialpr[j]
          end

          # initial guess:
          curr_x = initialpr

          # Probability of initial guess according to prior: π(Θ)
          curr_prior = logpdf(Distributions.MvNormal(curr_x, prior_σ), curr_x) # 10 allocations.
          if curr_prior == 0.0
            return "Probability of Prior too low"
          end
          if debug
            zeroparams = Array{Float64,1}()
            push!(zeroparams, curr_prior)
          end

          # Value of objective function at initial guess: Ln(Θ)
          curr_vals = objfun(curr_x) # 14 allocations / 9 kb

          # Probability of initial guess under proposal/Initialize a proposal probability. q(Θ'|Θ)
          curr_proposal_prob = 1

          while curr_it <= max_iterations
            for i =1:param_dim
              # Proposed next value, Θ'[i] - this is just *one* element.
              # Randomly generated conditional on current value, according to the proposal dist q(Θ'|Θ)
              next_x = curr_x #TODO: is this accepting every proposal?  Maybe.
              proposed = rand(proposal)[i]
              next_x[i] += proposed # This proposal is conditional - added noise to the current value.

              # Probability of proposed new value  Θ' under the prior: π(Θ')
              # This is *not* conditional on the current location
              # Using Log of normal PDF to avoid underflow.
              # Compute the prob of the entire new proposal
              next_prior = logpdf(Distributions.MvNormal(prior_μ, prior_σ), next_x) # 10 allocations / 1 kb

              # Probability of new value under proposal distribution:
              # This one *is* conditional on the current location.
              # Using Log of normal PDF to avoid underflow.
              next_proposal_prob = logpdf(Distributions.MvNormal(curr_x, pro_σ), next_x)

              # Value of objective at Proposal
              next_vals = objfun(next_x) #14 allocations / 9 kb

              # Difference in value of objectives
              val_diff = next_vals - curr_vals

              if val_diff == Inf || val_diff == 0.0  # keep track of times when this is too large.
                if val_diff == Inf
                  overflowcount += 1
                else val_diff == 0.0
                  underflowcount += 1
                end
              end
              #TODO: Note - maximizing or minimizing?  I want to get that right.  But it's irrelevant to the
              # algorithm since I can use -L instead of L

              logrho = minimum([val_diff+next_proposal_prob+next_prior-curr_prior-curr_proposal_prob,0.0]) #add the proposal.
              if debug
                trace[counter, 1] = logrho
                trace[counter, 2] = val_diff
                trace[counter, 3] = next_prior
                trace[counter, 4] = next_prior
                trace[counter, 5] = next_vals
                trace[counter, 6+i] = proposed
                for k =1:param_dim
                  allvals[counter,k] = next_x[k]
                end
              end
              if logrho >= 0 || rand() < exp(logrho)
                # Accepted Proposal
                curr_x = next_x
                curr_prior = next_prior
                curr_vals = next_vals
                curr_proposal_prob = next_proposal_prob
                accepted += 1
                if debug
                  param_accept[i]+= 1
                end
              end
              counter += 1
            end # of iteration over state elements.

            # Record the current parameter values whether they changed or not.
            # This is OUTSIDE the loop over parameter elements.
           for el = 1:param_dim
             @inbounds path[curr_it, el] = curr_x[el] # 4 allocations.
           end
            curr_it += 1
          end
          if debug
            return  path, overflowcount, underflowcount, accepted, trace, param_accept, allvals
          else
            return path, overflowcount, underflowcount, accepted
          end
        end # of MetropolisHastings()
    const nsims = 5000 #_000
    sim_vals, overcounter, undercounter, accept, tr, param_accept, allvals = MetropolisHastings(1000*ones(params), nsims)

    # Results to return:
    println("Fraction Accepted ", accept/nsims)
    println("Count of Numerical Overflow ", overcounter)
    println("Count of Underflow ", undercounter)

# remove latter two when not tracing results.
  return sim_vals, tr, param_accept, allvals

end # of LTE function


sims, tr, param_accept, allvals = LTE()

#=
p1 = plot(x=sims[:,1], Geom.histogram)
p2 = plot(x=sims[:,2], Geom.histogram)
p3 = plot(x=sims[:,3], Geom.histogram)
p4 = plot(x=sims[:,4], Geom.histogram)
p5 = plot(x=sims[:,5], Geom.histogram)
p6 = plot(x=sims[:,6], Geom.histogram)
p7 = plot(x=sims[:,7], Geom.histogram)

p62 = plot(x=sims[:,62], Geom.histogram)

=#
#probs = plot(x=zerop, Geom.histogram)




# for i = 1:size(sims, 2)
#   mode_n = mode(sims[:,i])
#   println("plot(x=sims[:,$i], Geom.histogram, Guide.xlabel(\" Parameter "*"$i"*" Mode $mode_n \"))")
# #  plot(x=sims[:,i], Geom.histogram)
# end





#
params_rec = zeros(size(sims,2))
for el in 1:size(sims, 2)
  tem = mean(sims[ convert(Int, size(sims,1)/4):end,el])
  println( tem)
  params_rec[el] = tem
end
#
# plot(x = collect(1:19), y=params_rec[1:19], Geom.line)

# function test()
#   x1 = rand(22, 130)
#   plot(x=mean(x1,1), Geom.histogram)
# end
# test()

#=
# for printing in the future:
for el in 1:size(Optim.minimizer(result3), 1)
  println(varcolnames[el+3], "  ", Optim.minimizer(result3)[el], " param: ", paramsymbs[el], " symbol number: ", el, "\n")
end
cnt = 0
for i=2:size(allvals,1)
  if allvals[i-1,1] != allvals[i, 1]
    println(allvals[i,1])
    cnt+=1
  end
end

=#
