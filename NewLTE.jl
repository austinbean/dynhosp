# Redoing the LTE.

using DataFrames
using Distributions
using StatsBase
using Plots  
plotlyjs()   # call PlotljJS backend to Plots.

#dat = readtable("/Users/austinbean/Google Drive/Simulation Results/dynhospsimresults.csv");
# This is just for testing purposes.  Not a real set of results.
#dat = readtable("/Users/austinbean/Desktop/dynhospsimulationresults2.csv");
dat = readtable("/Users/austinbean/Desktop/dynhospsimulationresults 10 03 2016 717pm.csv");

dat = convert(Array{Float64,2}, dat)

# TODO - right now this is not getting the constant revenue which is added by the medicaid patients.  

eq_opt = dat[:, 2:41];
neq_opt = dat[:,42:end];


function objfun(x::Vector; scale_fact = 1/10, inp1::Array{Float64,2}=scale_fact*eq_opt, inp2::Array{Float64,2}=scale_fact*neq_opt, diffmat::Array{Float64,2}=inp1-inp2)
  sum(min(diffmat*x, 0).^2)
end

function MetropolisHastings(initialpr::Vector,
                            max_iterations::Int64;
                            param_dim = length(initialpr),
                            pro_μ = zeros(param_dim),
                            pro_σ_scale::Float64 = 100.0,
                            pro_σ = pro_σ_scale*eye(param_dim),
                            proposal = Distributions.MvNormal(pro_μ, pro_σ),
                            prior_μ_scale::Float64 = 1000.0, # note that this should match the starting point, else the prob under the prior becomes 0
                            prior_μ = prior_μ_scale*ones(param_dim),
                            prior_σ_scale::Float64 = 100.0,
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
    tr= zeros(max_iterations*param_dim, 6+param_dim)
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

        proposed =  zeros(param_dim)
        disturb = rand(proposal)[i] #perturbation value.
        proposed[i] += disturb

        # Probability of proposed new value  Θ' under the *prior*: π(Θ')
        # This is *not* conditional on the current location
        # Using Log of normal PDF to avoid underflow.
        # Compute the prob of the entire new proposal under the prior π()
        next_prior = logpdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x + proposed) # could be done as logpdf(prior, next_x)

        # Probability of new value under the *proposal* q(x|y) distribution:
        # This one *is* conditional on the current location.
        # Using Log of normal PDF to avoid underflow.
        next_proposal_prob = logpdf(Distributions.MvNormal(curr_x, pro_σ), curr_x + proposed) # computing this isn't really necessary since it's not in the Hastings ratio.

        # Value of objective at Proposal
        next_vals = objfun(curr_x + proposed) #14 allocations / 9 kb

        # Difference in value of objectives
        val_diff = next_vals - curr_vals

        if val_diff == Inf || val_diff == 0.0  # keep track of times when this is too large.
          if val_diff == Inf
            overflowcount += 1 # basically can't have overflow since not taking exp() of objective
          else val_diff == 0.0
            underflowcount += 1 # this tracks when next_vals = current_vals, not underflow.  No exp() of objective anymore.
          end
        end
        # The probability under the proposal is symmetric - q(next_x|curr_x) = q(curr_x|next_x) - I can drop it.
        # Next line is log of Hastings ratio.
        logrho = minimum([val_diff+next_prior-curr_prior,0.0])
        if debug
          tr[counter, 1] = logrho
          tr[counter, 2] = val_diff
          tr[counter, 3] = next_prior
          tr[counter, 4] = next_proposal_prob # not using this for anything in the transition, but keeping track.
          tr[counter, 5] = next_vals
          tr[counter, 6+i] = disturb # trace records the disturbance to the value
          for k =1:param_dim
            @inbounds allvals[counter,k] = (curr_x + proposed)[k]
          end
        end
        if (logrho >= 0 || rand() < exp(logrho) )
          # Accepted Proposal
          curr_x[i] += disturb
          curr_prior = next_prior
          curr_vals = next_vals
          curr_proposal_prob = next_proposal_prob
          accepted += 1 # this is counting acceptances over all parameters.
          if debug
            param_accept[i] += 1
          end
          # Keep track of this too - should change the scaling factor.
          if logrho < 0
            if exp(logrho) == Inf
              overflowcount += 1
            end
          end
        end
        counter += 1
    end # of iteration over state elements.

      # Record the current parameter values whether they changed or not.
      # This is OUTSIDE the loop over parameter elements.
     for el = 1:param_dim
       # When proposal rejected, this should write out the old value.
       @inbounds path[curr_it, el] = curr_x[el] # 4 allocations.
     end
     curr_it += 1
  end
  if debug
    return  path, overflowcount, underflowcount, accepted, tr, param_accept, allvals
  else
    return path, overflowcount, underflowcount, accepted
  end
end # of MetropolisHastings()
const nsims = 1000 #_000
const paramsize = 40
sim_vals, overcounter, undercounter, accept, tr, param_accept, allvals = MetropolisHastings(1000*ones(paramsize), nsims)

# Results to return:
println("Fraction Accepted ", accept/nsims)
println("Count of Numerical Overflow ", overcounter)
println("Count of Underflow ", undercounter)






#=
# These columns just record what is where.
arr[1] = (alph1 = (beta^T)*outprob*dot(wtp_out[1:7], private[1:7])  )           # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 1
arr[2] = (alph2 = (beta^T)*outprob*dot(wtp_out[8:14], private[8:14])   )        # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 2
arr[3] = (alph3 = (beta^T)*outprob*dot(wtp_out[15:21], private[15:21]) )        # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 3
arr[4] = (gamma_1_385 = (beta^T)*outprob*(private[1]+medicaid[1])   )           # The next lines are patients summed over types.  Costs are treated as the same over Medicaid and privately insured.
arr[5] = (gamma_2_385 = (beta^T)*outprob*(private[8]+medicaid[8]) )             # I don't use any of the named terms - they just keep track.
arr[6] = (gamma_3_385 = (beta^T)*outprob*(private[15]+medicaid[15]) )
arr[7] = (gamma_1_386 = (beta^T)*outprob*(private[2]+medicaid[2]) )
arr[8] = (gamma_2_386 = (beta^T)*outprob*(private[9]+medicaid[9]) )
arr[9] = (gamma_3_386 = (beta^T)*outprob*(private[16]+medicaid[16]) )
arr[10] = (gamma_1_387 = (beta^T)*outprob*(private[3]+medicaid[3]) )
arr[11] = (gamma_2_387 = (beta^T)*outprob*(private[10]+medicaid[10]) )
arr[12] = (gamma_3_387 = (beta^T)*outprob*(private[17]+medicaid[17]) )
arr[13] = (gamma_1_388 = (beta^T)*outprob*(private[4]+medicaid[4]) )
arr[14] = (gamma_2_388 = (beta^T)*outprob*(private[11]+medicaid[11]) )
arr[15] = (gamma_3_388 = (beta^T)*outprob*(private[18]+medicaid[18]) )
arr[16] = (gamma_1_389 = (beta^T)*outprob*(private[5]+medicaid[5]) )
arr[17] = (gamma_2_389 = (beta^T)*outprob*(private[12]+medicaid[12]) )
arr[18] = (gamma_3_389 = (beta^T)*outprob*(private[19]+medicaid[19]) )
arr[19] = (gamma_1_390 = (beta^T)*outprob*(private[6]+medicaid[6]) )
arr[20] = (gamma_2_390 = (beta^T)*outprob*(private[13]+medicaid[13]) )
arr[21] = (gamma_3_390 = (beta^T)*outprob*(private[20]+medicaid[20]) )
arr[22] = (gamma_1_391 = (beta^T)*outprob*(private[7]+medicaid[7]) )
arr[23] = (gamma_2_391 = (beta^T)*outprob*(private[14]+medicaid[14]) )
arr[24] = (gamma_3_391 = (beta^T)*outprob*(private[21]+medicaid[21]) )
arr[25] = (patients385 = (beta^T)*outprob*(medicaid[1]+medicaid[8]+medicaid[15]))    # Count of patients at DRG 385
arr[26] = (patients386 = (beta^T)*outprob*(medicaid[2]+medicaid[9]+medicaid[16]))
arr[27] = (patients387 = (beta^T)*outprob*(medicaid[3]+medicaid[10]+medicaid[17]))
arr[28] = (patients388 = (beta^T)*outprob*(medicaid[4]+medicaid[11]+medicaid[18]))
arr[29] = (patients389 = (beta^T)*outprob*(medicaid[5]+medicaid[12]+medicaid[19]))
arr[30] = (patients390 = (beta^T)*outprob*(medicaid[6]+medicaid[13]+medicaid[20]))
arr[31] = (patients391 = (beta^T)*outprob*(medicaid[7]+medicaid[14]+medicaid[21]))
arr[32:end] = (alltrans = (beta^T)*outprob*transitions)

=#
