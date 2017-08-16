

"""
`MetropolisHastings`

Runs the LTE estimate on objfun above.
"""
function MetropolisHastings(initialpr::Vector,
                            max_iterations::Int64;
                            param_dim = length(initialpr),
                            pro_μ = zeros(param_dim),
                            pro_σ_scale::Float64 = 100.0,
                            pro_σ = pro_σ_scale*eye(param_dim),
                            proposal = Distributions.MvNormal(pro_μ, pro_σ),
                            prior_μ::Array{Float64,1} = [10.0, 1.0, 1.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412,  2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ],
                            prior_σ_scale::Float64 = 100.0,
                            prior_σ = prior_σ_scale*eye(param_dim),
                            prior = Distributions.MvNormal(prior_μ, prior_σ),
                            debug::Bool = false)
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
    if debug
      allvals[1,j] = initialpr[j]
    end
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
  preall = zeros(Float64, param_dim)
  proposed =  zeros(param_dim)
  while curr_it <= max_iterations
    for i =1:param_dim
        # Proposed next value, Θ'[i] - this is just *one* element.
        # Randomly generated conditional on current value, according to the proposal dist q(Θ'|Θ)

        ArrayZero(proposed)       # Clean up proposal 
        ArrayZero(preall)         # Clean up preallocated 
        rand!(proposal, preall)   # perturbation value. 
        proposed[i] += preall[i] # Changes just the one location 
        if !(find(proposed)==[i])
          println(proposed)
        end 
        if curr_it > 2 
          println("outside update: ", curr_prior, " iter ", curr_it)
        end 

        # Probability of proposed new value  Θ' under the *prior*: π(Θ')
        # This is *not* conditional on the current location
        # Using Log of normal PDF to avoid underflow.
        # Compute the prob of the entire new proposal under the prior π()
        pr1 = logpdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x)  # TODO - remove this line.  
        next_prior = logpdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x + proposed) # could be done as logpdf(prior, next_x)
        if pr1 == next_prior # this never happens.  
          println("priors equal? ", pr1, "  ", next_prior)
        end 
        if curr_it > 2
          print(" current pdf ", pdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x), " ")
          print(" proposed logpdf ", logpdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x + proposed), "  ")
          println("current pdf ", pdf(Distributions.MvNormal(prior_μ, prior_σ), curr_x + proposed))
        end 


        # Probability of new value under the *proposal* q(x|y) distribution:
        # This one *is* conditional on the current location.
        # Using Log of normal PDF to avoid underflow.
        # TODO - what is this line for?  
        next_proposal_prob = logpdf(Distributions.MvNormal(curr_x, pro_σ), curr_x + proposed) # computing this isn't really necessary since it's not in the Hastings ratio.

        # Value of objective at Proposal
        # This should be a scalar...
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
        # seems like proposals will almost only be accepted at random.
        # TODO - this is weird... next_prior and curr_prior are going to be small, no?  
        logrho = minimum([val_diff+next_prior-curr_prior,0.0])
        if i != 1
          if logrho > 0
            println("accepted!")
          end 
        end 
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
        # TODO - below this enables tracking acceptance 
        b1::Bool = false 
        b2::Bool = false 
        if logrho>=0
          b1 = true 
          println("accepted logrho ", logrho, " ", curr_it)
        end   
        if rand() < exp(logrho)
          b2 = true 
          println("accepted rand ", curr_it)
        end 
        # TODO - above this enables tracking acceptance.  
        if (b1)||(b2)
          # Accepted Proposal
          curr_x[i] += proposed[i]
          curr_prior = next_prior # this is a float.
          curr_vals = next_vals # this is a float?
          curr_proposal_prob = next_proposal_prob # this is a float.
          accepted += 1 # this is counting acceptances over all parameters.
          if debug
            param_accept[i] += 1
          end
          # Keep track of this too - should change the scaling factor.
          if logrho < 0
            # TODO - how can this condition ever be satisfied?  logrho < 0 but exp(logrho) = Inf?  
            if exp(logrho) == Inf
              println("some bullshit happpened")
              overflowcount += 1
            end
          end
        end
        if debug
          counter += 1
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
  # tr - a trace for debugging, including log of the function, the difference, the probability of the proposed val under the prior, the values, etc.
  # param_accept - Fraction of acceptances per parameter.
  # allvals - records the complete state at every round, so only one variable should change at any time.
  if debug
    return  path, overflowcount, underflowcount, accepted, tr, param_accept, allvals
  else
    return path, overflowcount, underflowcount, accepted
  end
end # of MetropolisHastings()