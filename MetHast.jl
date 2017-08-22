

"""
`MetropolisHastings`

Runs the LTE estimate on objfun above.

I suspect current issues are related to the probability only.  It must be that...
guess = [1000.0, 100.0, 10.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412, 2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ];
sim_vals, overcounter, undercounter, accepted = MetropolisHastings(guess, nsims, testfun) # no debugging output.

function tester()
  srand(1) # seed the generator to reproduce.
  nsims = 100_000
  guess = [1000.0, 100.0, 10.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412, 2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ];
  sim_vals, overcounter, undercounter, accepted = MetropolisHastings(guess, nsims, testfun, interimeq_opt, interimneq_opt, eq_const, neq_const; debug = false) # no debugging output.
  ans1 = GetModes(sim_vals);
  guess = [1000.0, 100.0, 10.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412, 2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ];
  ResultsPrint(ans1, guess)
end 
tester()

"""
function MetropolisHastings(initialpr::Vector,
                            max_iterations::Int64,
                            objfun::Function,  # take the function as an argument directly.
                            inp1::Array{Float64,2},   # arguments to objfun
                            inp2::Array{Float64,2},   # arguments to objfun
                            cons1::Array{Float64,2},  # arguments to objfun
                            cons2::Array{Float64,2};  # arguments to objfun
                            param_dim = length(initialpr),
                            pro_μ = zeros(param_dim),
                            pro_σ_scale::Float64 = 10.0,
                            pro_σ = pro_σ_scale*eye(param_dim),
                            proposal = Distributions.MvNormal(pro_μ, pro_σ),
                            prior_μ::Array{Float64,1} = [1000.0, 100.0, 10.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412,  2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ],
                            prior_σ_scale::Float64 = 10.0,
                            prior_σ = prior_σ_scale*eye(param_dim),
                            prior = Distributions.MvNormal(prior_μ, prior_σ))
  # Basics
  curr_it = 2
  overflowcount = 0
  underflowcount = 0
  priorzerocount = 0
  accepted = 0


  # Storing the values:
  path = zeros(max_iterations, param_dim) # 10 allocations / 50 MB
  for j = 1:param_dim
    path[1,j] = initialpr[j] #record initial guess
  end

  # initial guess:
  curr_x = initialpr

  # Probability of initial guess according to prior: π(Θ)
  curr_prior = logpdf(Distributions.MvNormal(curr_x, prior_σ), curr_x) # 10 allocations.
  if curr_prior == 0.0
    return "Probability of Prior too low"
  end


  # Value of objective function at initial guess: Ln(Θ)
  curr_vals = objfun(curr_x, inp1, inp2, cons1, cons2) # 14 allocations / 9 kb

  # Probability of initial guess under proposal/Initialize a proposal probability. q(Θ'|Θ)
  curr_proposal_prob = 1
  preall = zeros(Float64, param_dim)
  proposed =  zeros(param_dim)
  sign_check = zeros(Int64,param_dim)
  while curr_it <= max_iterations
    for i =1:param_dim
        # Proposed next value, Θ'[i] - this is just *one* element.
        # Randomly generated conditional on current value, according to the proposal dist q(Θ'|Θ)
        ArrayZero(proposed)       # Clean up proposal 
        ArrayZero(preall)         # Clean up preallocated 
        rand!(proposal, preall)   # perturbation value. 
        proposed[i] += preall[i] # Changes just the one location 
        if !(find(proposed)==[i])
          println(proposed) # this vector should be identically zero except at coordinate i.  
        end 

        # Probability of proposed new value  Θ' under the *prior*: π(Θ')
        # This is *not* conditional on the current location
        # Using Log of normal PDF to avoid underflow.
        # Compute the prob of the entire new proposal under the prior π()
        pr1 = logpdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x)    
        next_prior = logpdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x + proposed) # could be done as logpdf(prior, next_x)
        if pr1 == next_prior # this never happens.  
          println("priors equal? ", pr1, "  ", next_prior)
        end 
        # Probability of new value under the *proposal* q(x|y) distribution:
        # This one *is* conditional on the current location.
        # Using Log of normal PDF to avoid underflow.
        next_proposal_prob = logpdf(Distributions.MvNormal(curr_x, pro_σ), curr_x + proposed) 

        # Value of objective at Proposal
        next_vals = objfun(curr_x + proposed, inp1, inp2, cons1, cons2) 
        # Difference in value of objectives
        val_diff = next_vals - curr_vals
        if val_diff == 0.0 # this almost never happens.  
          println("Value difference is 0.0 ", i)
        end 

        if val_diff == Inf || val_diff == 0.0  # keep track of times when this is too large.
          if val_diff == Inf
            overflowcount += 1 # basically can't have overflow since not taking exp() of objective
          else val_diff == 0.0
            underflowcount += 1 # this tracks when next_vals = current_vals, not underflow.  No exp() of objective anymore.
          end
        end
        # The probability under the proposal is symmetric - q(next_x|curr_x) = q(curr_x|next_x) - I can drop it.
        # NOTE - the proposal distribution is symmetric so the following quantities are identical. 
        prop_cond_curr = logpdf(Distributions.MvNormal(curr_x + proposed, pro_σ), curr_x) # this one subtracted
        curr_cond_prop = logpdf(Distributions.MvNormal(curr_x, pro_σ), curr_x + proposed) # this one added
        acpt::Bool = false
        p::Float64 = min(exp(next_vals-curr_vals)*(next_prior/curr_prior)*(curr_cond_prop/prop_cond_curr),1)
        acpt = StatsBase.sample([true, false], StatsBase.weights([p, 1-p]))
        # The following rejects proposals outside the bounds of the parameter space.
        if (i ==28)||(i == 31)||(i ==32)   # params that must be less than 0
          if (curr_x[i]+proposed[i]>0)   # check for parameter greater than zero and reject if so 
            acpt = false                   # reject anything greater than 0
            sign_check[i] += 1
          end 
        else                               # params which cannot be less than 0
          if (curr_x[i]+proposed[i] < 0)   # check for parameter less than zero and reject if so 
            acpt = false                   # reject anything less than 0
            sign_check[i] += 1
          end
        end 
        if acpt
          # Accepted Proposal
          curr_x[i] += proposed[i]
          curr_prior = next_prior # this is a scalar
          curr_vals = next_vals # this is a float, equal to the objective function at the present parameters.
          curr_proposal_prob = next_proposal_prob # this is a float.
          accepted += 1 # this is counting acceptances over all parameters.
        end
    end # of iteration over state elements.
    # Record the current parameter values whether they changed or not.
    # This is OUTSIDE the loop over parameter elements.
    for el = 1:param_dim
      @inbounds path[curr_it, el] = curr_x[el] # When proposal rejected, this should write out the old value.
    end
    curr_it += 1
  end
  # Path - the accepted values from one complete round of proposals to all parameters
  # Overflowcount - number of times the function overflowed
  # Underflowcount - number of times the function underflowed
  # Accepted - Counts the number of times the proposal was accepted (one parameter at a time)
  for k = 1:size(sign_check,1)
    println(sign_check[k])
  end 
  return path, overflowcount, underflowcount, accepted
end # of MetropolisHastings()
