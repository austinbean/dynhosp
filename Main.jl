# Main

# Packages
using DataFrames
using DataArrays
using Distributions

# Include necessary functions
include("/Users/austinbean/Desktop/dynhosp/combgen.jl")
include("/Users/austinbean/Desktop/dynhosp/nckr.jl")
include("/Users/austinbean/Desktop/dynhosp/probfinder.jl")
include("/Users/austinbean/Desktop/dynhosp/probfind2.jl")
include("/Users/austinbean/Desktop/dynhosp/tuplefinder.jl")


# Import Data
dataf = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
notmissing = findin(isna(dataf[:fipscode]), false);
dataf = dataf[notmissing, :];

# regcoeffs = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/ProbabilityCoeffs.csv", header = true);

simdata = readtable("/Users/austinbean/Desktop/dynhosp/Simulated Choice Probs.csv", header = true);
sim_f = DataFrame(simdata)



### Collect Basic Information ###

# locates starting and ending points of all separate fipscode indices
allindices = [ (x, findfirst(dataf[:fipscode], x), findlast(dataf[:fipscode], x)) for x in unique(dataf[:fipscode]) ];
# Next one also works in case I decide an array of arrays is better than an array of tuples
#indices = [ [x, findfirst(df[:fipscode], x), findlast(df[:fipscode], x)] for x in unique(df[:fipscode]) ]
# Next one stores all the years appearing in each fips code
yearins = [ [x; findfirst(dataf[:fipscode], x); findlast(dataf[:fipscode], x ); unique( dataf[findfirst(dataf[:fipscode], x):findlast(dataf[:fipscode], x ) , :year]  ) ] for x in unique(dataf[:fipscode])  ]


######  Simulation Section Starts Here #######

β = 0.95;
T = 100;
choices = 4;
α₂ = 0.07;  # Fraction of patients admitted to NICU lev 2 on average (PA Data)
α₃ = 0.13; # Fraction of patients admitted to NICU Lev 3 on average (PA Data)

# Parameters of the shock distribution
# This is wrong - I don't draw shocks from this, I draw conditional means,
# conditional on being the maximum.
dist_μ = 0;
dist_σ = 1;
dist_ξ = 0;
srand(123)
d = GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# I do need the constant:
γ = eulergamma;

# Fix the simulation at Level 1, just to start:
# but then what do we do when he switches?
# Record the firm's own state as the tuple ownstate with (0,0) - lev 1 (1, 0) - lev 2 (0,1) - lev 3
ownstate = (1,0)
total_hosp = 8;
level1 = 3;
level2 = 2;
level3 = 3;
# Three consecutive columns in configs storing future market arrangements
future_config1 = 4;
future_config2 = 5;
future_config3 = 6;

# Possible choices:
choices_1 = [(1,1), (1,2), (1,3), (1,4)];
choices_2 = [(2,1), (2,2), (2,3), (2,4)];
choices_3 = [(3,1), (3,2), (3,3), (3,4)];



# Actual entry probabilities will be substituted in later.
entryprobs = [0.99, 0.004, 0.001, 0.005] # [No entry, level1, level2, level3] - not taken from anything, just imposed.
entrants = [0, 1, 2, 3]

# Note: since I know the simulation periods in advance, I do *not* need to resize dynamically, just fill in.
#state_history = [ownstate 1 1 level1 level2 level3]

#=
Next thing - 04 05 16
• Must keep track of who is where.
• Add the whole set of hospital neighbors to each year-record
• How far?  Clearly 50 miles is too far - that will be up to 52 neighbors: far too many
• Continue with this part of the simulation without worrying about the error.  that can be fixed
• Need to write a function to get the logit probs, now that it depends on distance.

• Draw uniform probs instead of shocks...  (???)
• Choose the action based on those probs
• Perturb policies as already imagined: change by ϵ the probs of actions.

- Also: need to finally write the logit probability function for unobserved state
combinations to get those probabilities.
- Or does it make sense to do something like a local linear regression over the
unobserved state space elements?
=#


start = 2;

for y in 1:size(yearins)[1]
	market_start = yearins[y][2]
	market_end = yearins[y][3]
	market_frame = dataf[market_start:market_end, :]
	for year in yearins[y][4:end]
		year_frame = market_frame[(market_frame[:, :year].==year), :]
		level1 = year_frame[1,:level1_hospitals0]
		level2 = year_frame[1,:level2solo_hospitals0]
		level3 = year_frame[1,:level3_hospitals0]
		fids = unique(year_frame[:fid])
		all_hosp_probs = zeros(size(fids)[1], 11)
		for fid in 1:size(fids)[1]
			el = fids[fid]
			a = year_frame[:fid].== el
			all_hosp_probs[fid, 1:end] =hcat(el, year_frame[a, :act_int], year_frame[a, :act_solo], year_frame[a, :choicenum0], year_frame[a, :pr_ch_0], year_frame[a, :choicenum1], year_frame[a, :pr_ch_1], year_frame[a, :choicenum2], year_frame[a, :pr_ch_2], year_frame[a, :choicenum3], year_frame[a, :pr_ch_3] )
		end
		state_history = [ownstate 1 1 level1 level2 level3; zeros(T, 6)]
		for i = start:T+1
			rand_action = zeros(size(fids)[1], 3)
			for hosp in 1:size(all_hosp_probs)[1]
				if ((all_hosp_probs[hosp,2] == 0) & (all_hosp_probs[hosp,3] == 0)) # level 1
					pairs = hcat(transpose([all_hosp_probs[hosp,4] all_hosp_probs[hosp,6] all_hosp_probs[hosp,8] all_hosp_probs[hosp,10]]), choices_1, transpose([all_hosp_probs[hosp,5] all_hosp_probs[hosp,7] all_hosp_probs[hosp,9] all_hosp_probs[hosp,11]])  )
				elseif ((all_hosp_probs[hosp,2] == 1) & (all_hosp_probs[hosp,3] == 0))
					pairs = hcat(transpose([all_hosp_probs[hosp,4] all_hosp_probs[hosp,6] all_hosp_probs[hosp,8] all_hosp_probs[hosp,10]]), choices_2, transpose([all_hosp_probs[hosp,5] all_hosp_probs[hosp,7] all_hosp_probs[hosp,9] all_hosp_probs[hosp,11]])  )
				elseif ((all_hosp_probs[hosp,2] == 1) & (all_hosp_probs[hosp,3] == 0))
					pairs = hcat(transpose([all_hosp_probs[hosp,4] all_hosp_probs[hosp,6] all_hosp_probs[hosp,8] all_hosp_probs[hosp,10]]), choices_3, transpose([all_hosp_probs[hosp,5] all_hosp_probs[hosp,7] all_hosp_probs[hosp,9] all_hosp_probs[hosp,11]])  )
				end
				weights = WeightVec(Array{Float64}(pairs[:, 3]))
				draws = sample( pairs[:,2], weights)
				poutcomes = [draws, pairs[findfirst(pairs[:,2], draws), 3]]
				rand_action[hosp,:] = hcat(all_hosp_probs[hosp, 1], poutcomes[1], poutcomes[2] )
			end
			# Entry draws
			entrypairs = hcat(entrants, entryprobs)
			entrantsp = WeightVec(entryprobs)
			newentrant = sample(entrants, entrantsp)
			entrantout = [newentrant, entrypairs[findfirst(entrypairs[:,1], newentrant), 2]]'
			# Transitions:
			total_hosp = level1 + level2 + level3
			state_history[i, :] =  [ownstate probstate probstate*state_history[i-1,3] level1 level2 level3]
				#state_history = vcat(state_history, [ownstate probstate probstate*state_history[end,3] level1 level2 level3])
# Now we need to track the aggregate state.




			elseif ownstate == (1,0)
				if (!((level1 == prev1) & (level2 == prev2) & (level3 == prev3))) # Won't just be relevant at *start*
					global pvect
					if findprob(dataf, total_hosp, level1, level2, level3) != ProbException
						pvect = findprob(dataf, total_hosp, level1, level2, level3)
					elseif findprob(dataf, total_hosp, level1, level2, level3) == ProbException
						if findprob(simdata, total_hosp, level1, level2, level3) != ProbException
							pvect = findprob(simdata, total_hosp, level1, level2, level3)
						else
							println("Probability Exception at Configuration (L1, L2, L3) = ", level1, level2, level3)
							println("Market configuration not found in actual or simulated data")
							break
						end
					end
				end
				# shockdraw = γ -log(pvect_2[5:8]).*pvect_2[5:8] # This is wrong - correct later
				choices_2 = pvect[5:8]
				configs =sortrows( nckrexen([max(level1, 0) max(level2-1, 0) max(level3, 0)]), by = x -> (x[4], x[5], x[6]) ) #future market configs ignoring the level 1 hospital
				futures = tuplefinder(configs, future_config1, future_config2, future_config3) # potential arrangements of the market, minus the current hospital
				transitions =hcat( [configs[:,4] configs[:,5] configs[:,6]], prod(broadcast(^,hcat(pvect', entryprobs'), configs[:, 7:22] ), 2))
				stateholder = hcat(futures, zeros(size(futures)[1], 1))
				for k = 1:size(transitions)[1]
					temptup = (transitions[k,1], transitions[k,2], transitions[k,3])
					stateholder[findfirst(stateholder[:,1], temptup), 2] += transitions[k, 4]
				end
				stweights =WeightVec(Array{Float64}(stateholder[:,2]))
				next_state = sample(stateholder[:,1], stweights)
				probstate = stateholder[findfirst(stateholder[:,1], next_state), 2]
				lp21 = log(choices_2[1]) + shockdraw[1]
				lp22 = log(choices_2[2]) + shockdraw[2]
				lp23 = log(choices_2[3]) + shockdraw[3]
				lp2EX = log(choices_2[4]) + shockdraw[4]
				max_choice = indmax([lp21, lp22, lp23, lp2EX])
				if max_choice == 1
					ownstate = (0,0)
					prev1 = level1
					prev2 = level2
					prev3 = level3
					level1 = next_state[1] + 1
					level2 = next_state[2]
					level3 = next_state[3]
				elseif max_choice == 2
					ownstate = (1,0)
					prev1 = level1
					prev2 = level2
					prev3 = level3
					level1 = next_state[1]
					level2 = next_state[2] + 1
					level3 = next_state[3]
				elseif max_choice == 3
					ownstate = (0,1)
					prev1 = level1
					prev2 = level2
					prev3 = level3
					level1 = next_state[1]
					level2 = next_state[2]
					level3 = next_state[3] + 1
				elseif max_choice == 4
					ownstate = ("EX", "EX")
					level1 = next_state[1]
					level2 = next_state[2]
					level3 = next_state[3]
				else
					println("F-ed Up")
				end
				total_hosp = level1 + level2 + level3
				state_history[i, :] =  [ownstate probstate probstate*state_history[i-1,3] level1 level2 level3]

		#	state_history = vcat(state_history, [ownstate probstate probstate*state_history[end,3] level1 level2 level3])
			elseif ownstate == (0,1)
				if (!((level1 == prev1) & (level2 == prev2) &(level3 == prev3)) )
					global pvect
					if findprob(dataf, total_hosp, level1, level2, level3) != ProbException
						pvect = findprob(dataf, total_hosp, level1, level2, level3)
					elseif findprob(dataf, total_hosp, level1, level2, level3) == ProbException
						if findprob(simdata, total_hosp, level1, level2, level3) != ProbException
							pvect = findprob(simdata, total_hosp, level1, level2, level3)
						else
							println("Probability Exception at Configuration (L1, L2, L3) = ", level1, level2, level3)
							println("Market configuration not found in actual or simulated data")
							break
						end
					end
				end
				# shockdraw = γ -log(pvect_3[9:12]).*pvect_3[9:12] # this is wrong - correct later.
				choices_3 = pvect[9:12]
				configs =sortrows( nckrexen([max(level1, 0) max(level2, 0) max(level3-1, 0)]), by = x -> (x[4], x[5], x[6]) ) #future market configs ignoring the level 1 hospital
				futures = tuplefinder(configs, future_config1, future_config2, future_config3) # potential arrangements of the market, minus the current hospital
				transitions =hcat( [configs[:,4] configs[:,5] configs[:,6]], prod(broadcast(^,hcat(pvect', entryprobs'), configs[:, 7:22] ), 2))
				stateholder = hcat(futures, zeros(size(futures)[1], 1))
				for k = 1:size(transitions)[1]
					temptup = (transitions[k,1], transitions[k,2], transitions[k,3])
					stateholder[findfirst(stateholder[:,1], temptup), 2] += transitions[k, 4]
				end
				stweights =WeightVec(Array{Float64}(stateholder[:,2]))
				next_state = sample(stateholder[:,1], stweights)
				probstate = stateholder[findfirst(stateholder[:,1], next_state), 2]
				lp31 = log(choices_3[1]) + shockdraw[1]
				lp32 = log(choices_3[2]) + shockdraw[2]
				lp33 = log(choices_3[3]) + shockdraw[3]
				lp3EX = log(choices_3[4]) + shockdraw[4]
				max_choice = indmax([lp31, lp32, lp33, lp3EX])
				if max_choice == 1
					ownstate = (0,0)
					prev1 = level1
					prev2 = level2
					prev3 = level3
					level1 = next_state[1] + 1
					level2 = next_state[2]
					level3 = next_state[3]
				elseif max_choice == 2
					ownstate = (1,0)
					prev1 = level1
					prev2 = level2
					prev3 = level3
					level1 = next_state[1]
					level2 = next_state[2] + 1
					level3 = next_state[3]
				elseif max_choice == 3
					ownstate = (0,1)
					prev1 = level1
					prev2 = level2
					prev3 = level3
					level1 = next_state[1]
					level2 = next_state[2]
					level3 = next_state[3] + 1
				elseif max_choice == 4
					ownstate = ("EX", "EX")
					level1 = next_state[1]
					level2 = next_state[2]
					level3 = next_state[3]
				else
					println("F-ed Up")
				end
				total_hosp = level1 + level2 + level3
				state_history[i, :] =  [ownstate probstate probstate*state_history[i-1,3] level1 level2 level3]

		#		state_history = vcat(state_history, [ownstate probstate probstate*state_history[end,3] level1 level2 level3])
			elseif ownstate == ("EX", "EX")
				state_history[i, :] =  [ownstate 1 1 level1 level2 level3]

				state_history = vcat(state_history, [ownstate 1 1 level1 level2 level3])
			end
		end
	end
end

# Extract the count of visits to various states in this section
# Lev 1: (0,0)
# Lev 2: (1,0)
# Lev 3: (0,1)

function pairdiff(x::Tuple, y::Tuple)
	return x[1] - y[1], x[2] - y[2]
end



demand = [ 100, 150, 200]
previous_own = state_history[1,1]
previous_agg = state_history[1,2:4]

for row in 2:size(state_history)[1]
	discount = β^row
	println(state_history[row,:])
	# Record changes in the aggregate state and the own state
 	current_own = state_history[row, 1]
	current_agg = state_history[row, 2:4]
	change_own = pairdiff(current, previous)
	change_agg = current_agg - previous_agg # do I care?
	# Record visits to possible aggregate states

#=

- Use the histogram function.  See this: http://stackoverflow.com/questions/21172027/count-instances-of-each-unique-integer-in-a-vector-in-1-line-of-code

=#
	previous = state_history[row, 1] # records the current state as "previous for the next round"
end












###
