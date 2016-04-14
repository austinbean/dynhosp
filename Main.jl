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
include("/Users/austinbean/Desktop/dynhosp/LogitEst.jl")


# Import Data
dataf = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
notmissing = findin(isna(dataf[:fipscode]), false);
dataf = dataf[notmissing, :];

regcoeffs = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Choice Model.csv", header = true);

#simdata = readtable("/Users/austinbean/Desktop/dynhosp/Simulated Choice Probs.csv", header = true);
#sim_f = DataFrame(simdata)



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



type Hosp
	fid::Int
	name::AbstractString
	level::Tuple
	fips::Int
	lat::Float64
	long::Float64
	n05::Array
	n515::Array
	n1525::Array
	choices::Array
	probs::WeightVec
end

h1 = Hosp(dataf[1,:fid], dataf[1,:facility], (dataf[1,:act_int], dataf[1,:act_solo]), dataf[1,:fipscode], dataf[1, :v15], dataf[1, :v16],     [dataf[1,:lev105], dataf[1,:lev205], dataf[1,:lev305]], [dataf[1,:lev1515], dataf[1,:lev2515], dataf[1,:lev3515]], [dataf[1,:lev11525], dataf[1,:lev21525], dataf[1,:lev31525]], [1, 2, 3, 4], WeightVec([0.1, 0.1, 0.1, 0.1])     )

type Market
	fips::Int
	lev1::Int
	lev2::Int
	lev3::Int
	config::Array{Hosp}
end

# Action codes:
#=
1 "To 3 from 1"
2 "To 2 from 1"
3 "To 2 from 3"
4 "To 1 from 3"
5 "To 1 from 2"
6 "To 3 from 2"
7 "Enter at 1"
8 "Enter at 2"
9 "Enter at 3"
10 "Do Nothing"
11 "Exit"
=#

# Problem: hospitals may be closer than 25 miles but in different counties - they
# will show up as being in different markets.

start = 2;
neighbors_start = 108;
fields = 6;


for y in 1:size(yearins)[1]
	market_start = yearins[y][2]
	market_end = yearins[y][3]
	market_frame = dataf[market_start:market_end, :]
		for year in yearins[y][4:end]
			println(year)
			year_frame = market_frame[(market_frame[:, :year].==year), :]
			level1 = year_frame[1,:level1_hospitals0]
			level2 = year_frame[1,:level2solo_hospitals0]
			level3 = year_frame[1,:level3_hospitals0]
			fids = unique(year_frame[:fid])
			# What do I want to track over the whole history? fid, solo state, int state, action chosen, probability of choice, demand.  Aggregate: prob, levels.
			# Also think forward: demand realized.
			state_history = [zeros(1, fields*size(fids)[1]) 1 level1 level2 level3; zeros(T, fields*(size(fids)[1]) + 4)]
			for i = start:T+1
				println(i)
						for fid in 1:size(fids)[1] # this has to be handled separately for each hospital, due to the geography issue
							el = fids[fid] # the dataframe is mutable.
							a = (year_frame[:fid].==el)
							if  ((year_frame[a,:act_int][1], year_frame[a,:act_solo][1]) == (0,0)) # level 1, actions:
								probs1 = logitest((0,0), level1, level2, level3, convert(Array, [year_frame[a,:lev105]; year_frame[a,:lev205]; year_frame[a,:lev305]; year_frame[a,:lev1515]; year_frame[a,:lev2515]; year_frame[a,:lev3515]; year_frame[a,:lev11525]; year_frame[a,:lev21525]; year_frame[a,:lev31525]]) )
								# all_hosp_probs[fid] = hcat(el, year_frame[a, :act_int], year_frame[a,:act_solo], 10, probs[1], 2, probs[2], 1, probs[3], 11, probs[4])
								# Reassign action choices -
								year_frame[a, :choicenum0] = 10
								year_frame[a, :pr_ch_0] = probs1[1]
								year_frame[a, :choicenum1] = 2
								year_frame[a, :pr_ch_1] = probs1[2]
								year_frame[a, :choicenum2] = 1
								year_frame[a, :pr_ch_2] = probs1[3]
								year_frame[a, :choicenum3] = 11
								year_frame[a, :pr_ch_3] = probs1[4]
								# Draw action:
								action1 = sample([10, 2, 1, 11] ,WeightVec([probs1[1], probs1[2], probs1[3], probs1[4]]))
								# Change things to reflect the action chosen:
								# record the probability of the action taken.
								if action1 == 10
									chprob = probs1[1]
									# no change to state in the aggregate
									# no change to distance-state of others
								elseif action1 == 2
									chprob = probs1[2]
									year_frame[a,:act_int] = 1
									year_frame[a,:act_solo] = 0
								elseif action1 == 1
									chprob = probs1[3]
									year_frame[a,:act_int] = 0
									year_frame[a,:act_solo] = 1
								elseif action1 == 11
									chprob = probs1[4]
									year_frame[a,:act_int] = -999
									year_frame[a,:act_solo] = -999
								else
									println("Fail")
									println("Bad Action Chosen by", el)
									break
								end
								# Set own distance counts to 0 for all categories
								(year_frame[a,:lev105], year_frame[a,:lev205], year_frame[a,:lev305], year_frame[a,:lev1515], year_frame[a,:lev2515], year_frame[a,:lev3515], year_frame[a,:lev11525], year_frame[a,:lev21525], year_frame[a,:lev31525]) = zeros(1,9)
								# This is the logical place to recompute the demand.
								# write out state values - in blocks:
								state_history[i, (fid-1)*fields + 1] = el # Change
								state_history[i, (fid-1)*fields + 2] = year_frame[a,:act_int][1] # can't assign tuple to int
								state_history[i, (fid-1)*fields + 3] = year_frame[a,:act_solo][1] # can't assign tuple to int
								state_history[i, (fid-1)*fields + 4] = chprob
								state_history[i, (fid-1)*fields + 5] = action1
								#state_history[i, (fid-1)*fields + 5] = demand
							elseif ((year_frame[a,:act_int][1], year_frame[a,:act_solo][1]) == (1,0)) #level 2, actions:
								# Can't evaluation as nexti here if they are initialized to negative 1
								probs2 = logitest((1,0), level1, level2, level3, convert(Array, [year_frame[a,:lev105]; year_frame[a,:lev205]; year_frame[a,:lev305]; year_frame[a,:lev1515]; year_frame[a,:lev2515]; year_frame[a,:lev3515]; year_frame[a,:lev11525]; year_frame[a,:lev21525]; year_frame[a,:lev31525]]) )
								#all_hosp_probs[fid] = hcat(el, year_frame[a, :act_int], year_frame[a,:act_solo], 5, probs[1], 10, probs[2], 6, probs[3], 11, probs[4])
								year_frame[a, :choicenum0] = 5
								year_frame[a, :pr_ch_0] = probs2[1]
								year_frame[a, :choicenum1] = 10
								year_frame[a, :pr_ch_1] = probs2[2]
								year_frame[a, :choicenum2] = 6
								year_frame[a, :pr_ch_2] = probs2[3]
								year_frame[a, :choicenum3] = 11
								year_frame[a, :pr_ch_3] = probs2[4]
								# Action:
								action2 = sample([5, 10, 6, 11] ,WeightVec([probs2[1], probs2[2], probs2[3], probs2[4]]))
								if action2 == 5
									chprob = probs2[1]
									year_frame[a,:act_int] = 0
									year_frame[a,:act_solo] = 0
								elseif action2 == 10
									chprob = probs2[2]
									# no change to state
								elseif action2 == 6
									chprob = probs2[3]
									year_frame[a,:act_int] = 1
									year_frame[a,:act_solo] = 0
								elseif action2 == 11
									chprob = probs2[4]
									year_frame[a,:act_int] = -999
									year_frame[a,:act_solo] = -999
								else
									println("Fail")
									println("Bad Action Chosen by", el)
									break
								end
								# Set own distance counts to 0 for all categories
								(year_frame[a,:lev105], year_frame[a,:lev205], year_frame[a,:lev305], year_frame[a,:lev1515], year_frame[a,:lev2515], year_frame[a,:lev3515], year_frame[a,:lev11525], year_frame[a,:lev21525], year_frame[a,:lev31525]) = zeros(1,9)
								# write out state values - in blocks:
								state_history[i, (fid-1)*fields + 1] = el
								state_history[i, (fid-1)*fields + 2] = year_frame[a,:act_int][1]
								state_history[i, (fid-1)*fields + 3] = year_frame[a,:act_solo][1]
								state_history[i, (fid-1)*fields + 4] = chprob
								state_history[i, (fid-1)*fields + 5] = action2
								#state_history[i, (fid-1)*fields + 6] = demand
							elseif  ((year_frame[a,:act_int][1], year_frame[a,:act_solo][1]) == (0,1)) #level 3, actions:
								probs3 = logitest((0,1), level1, level2, level3, convert(Array, [year_frame[a,:lev105]; year_frame[a,:lev205]; year_frame[a,:lev305]; year_frame[a,:lev1515]; year_frame[a,:lev2515]; year_frame[a,:lev3515]; year_frame[a,:lev11525]; year_frame[a,:lev21525]; year_frame[a,:lev31525]]) )
								#all_hosp_probs[fid] = hcat(el, year_frame[a, :act_int], year_frame[a,:act_solo], 4, probs[1], 3, probs[2], 10, probs[3], 11, probs[4])
								year_frame[a, :choicenum0] = 4
								year_frame[a, :pr_ch_0] = probs3[1]
								year_frame[a, :choicenum1] = 3
								year_frame[a, :pr_ch_1] = probs3[2]
								year_frame[a, :choicenum2] = 10
								year_frame[a, :pr_ch_2] = probs3[3]
								year_frame[a, :choicenum3] = 11
								year_frame[a, :pr_ch_3] = probs3[4]
								# Action:
								action3 = sample([4, 3, 10, 11] ,WeightVec([probs3[1], probs3[2], probs3[3], probs3[4]]))
								if action3 == 3
									chprob = probs3[2]
									year_frame[a, :act_int] = 1
									year_frame[a, :act_solo] = 0
								elseif action3 == 4
									chprob = probs3[1]
									year_frame[a, :act_int] = 0
									year_frame[a, :act_solo] =0
								elseif action3 == 10
									chprob = probs3[3]
									# no change to state
								elseif action3 == 11
									chprob = probs3[4]
									year_frame[a, :act_int] = -999
									year_frame[a, :act_solo] = -999
								else
									println("Fail")
									println("Bad Action Chosen by", el)
									break
								end
								# Set own distance counts to 0 for all categories
								(year_frame[a,:lev105], year_frame[a,:lev205], year_frame[a,:lev305], year_frame[a,:lev1515], year_frame[a,:lev2515], year_frame[a,:lev3515], year_frame[a,:lev11525], year_frame[a,:lev21525], year_frame[a,:lev31525]) = zeros(1,9)
								# write out state values - in blocks:
								state_history[i, (fid-1)*fields + 1] = el
								state_history[i, (fid-1)*fields + 2] = year_frame[a,:act_int][1]
								state_history[i, (fid-1)*fields + 3] = year_frame[a,:act_solo][1]
								state_history[i, (fid-1)*fields + 4] = chprob
								state_history[i, (fid-1)*fields + 5] = action3
								#state_history[i, (fid-1)*fields + 6] = demand
							elseif ((year_frame[a,:act_int][1], year_frame[a,:act_solo][1]) == (-999,-999)) # has exited.
								# No new actions to compute, but record.
								state_history[i, (fid-1)*fields + 1] = el
								state_history[i, (fid-1)*fields + 2] = year_frame[a,:act_int][1]
								state_history[i, (fid-1)*fields + 3] = year_frame[a,:act_solo][1]
								state_history[i, (fid-1)*fields + 4] = 1 # exit is absorbing, so the choice prob is always 1
								state_history[i, (fid-1)*fields + 5] = 0 # no action is taken.
								# Set own distance counts to 0 for all categories
								(year_frame[a,:lev105], year_frame[a,:lev205], year_frame[a,:lev305], year_frame[a,:lev1515], year_frame[a,:lev2515], year_frame[a,:lev3515], year_frame[a,:lev11525], year_frame[a,:lev21525], year_frame[a,:lev31525]) = zeros(1,9)
							end
						end
						# Count facilities by distance and map results to neighboring hospitals
						# here the issue is that, for hospitals in neighboring counties, we haven't set the levX_YZ values to 0
						for i = 1:size(year_frame)[1]
							own_fac = (year_frame[i, :act_solo], year_frame[i, :act_int])
							for j = neighbors_start:(2):size(year_frame)[2]
								if (!isna(year_frame[i,j])) & (!isna(year_frame[i,j+1]))
									if (year_frame[i, j+1] >0) & (year_frame[i,j+1] < 5)
										if own_fac == (0,0)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev105] += 1
											println("updated facility")
										elseif own_fac == (1,0)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev205] += 1
											println("updated facility")
										elseif own_fac == (0,1)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev305] += 1
											println("updated facility")
										elseif own_fac == (-999,-999)
											# do nothing - firm exited.
										else
											println("Bad Own Facility Code", i, j, own_fac)
										end
									elseif (year_frame[i, j+1] >5) & (year_frame[i,j+1] < 15)
										if own_fac == (0,0)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev1515] += 1
											println("updated facility")
										elseif own_fac == (1,0)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev2515] += 1
											println("updated facility")
										elseif own_fac == (0,1)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev3515] += 1
											println("updated facility")
										elseif own_fac == (-999,-999)
											# do nothing - firm exited.
										else
											println("Bad Own Facility Code", i, j, own_fac)
										end
									elseif (year_frame[i,j+1] > 15) & (year_frame[i,j+1] < 25)
										if own_fac == (0,0)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev11525] += 1
										elseif own_fac == (1,0)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev21525] += 1
										elseif own_fac == (0,1)
											dataf[ (dataf[:fid].== year_frame[i,j])&(dataf[:year].== year), :lev31525] += 1
										elseif own_fac == (-999,-999)
											# do nothing - firm exited.
										else
											println("Bad Own Facility Code", i, j, own_fac)
										end
									else
										println("Bad Distance at frame ", i, j)
									end
								end
							end
						end
						# Sum the levels for next period:
						level1 = 0; level2 = 0; level3 = 0;
						for row in 1:size(year_frame)[1]
							if  (year_frame[row,:act_solo], year_frame[row,:act_int]) ==(0,0) # count all firms at each row in year_frame
								level1 += 1
							elseif (year_frame[row,:act_solo], year_frame[row,:act_int]) ==(1,0)
								level2 += 1
							elseif (year_frame[row,:act_solo], year_frame[row,:act_int]) ==(0,1)
								level3 += 1
							end
						end
						# Entry draw
						entrypairs = hcat(entrants, entryprobs)
						entrantsp = WeightVec(entryprobs)
						newentrant = sample(entrants, entrantsp)
						entrantout = [newentrant, entrypairs[findfirst(entrypairs[:,1], newentrant), 2]]'
						if newentrant> 0
							# Update fid count in here too.
							# Find a location for this guy too
							if newentrant == 1
								level1 += 1
							elseif newentrant == 2
								level2 += 1
							elseif newentrant == 3
								level3 += 1
						end
						# tracking the state history is fine for developing, but I really need to track every firm's history.
						total = level1 + level2 + level3

						# If someone entered here, I need to add a new fid.
						state_history[i, fid*fields+1] = level1 ;
						state_history[i, fid*fields+2] = level2 ;
						state_history[i, fid*fields+3] = level3 ;
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
