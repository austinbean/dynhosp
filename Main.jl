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
include("/Users/austinbean/Desktop/dynhosp/Distance.jl")


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



#=

• Choose the action based on those probs
• Perturb policies as already imagined: change by ϵ the probs of actions.


 Problem: hospitals may be closer than 25 miles but in different counties - they
will show up as being in different markets.
=#

# Data Types - currently unused.

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


######  Simulation Section Starts Here #######

β = 0.95;
T = 100;
choices = 4;
α₂ = 0.07;  # Fraction of patients admitted to NICU lev 2 on average (PA Data)
α₃ = 0.13; # Fraction of patients admitted to NICU Lev 3 on average (PA Data)
# Actual entry probabilities will be substituted in later.
entryprobs = [0.99, 0.004, 0.001, 0.005] # [No entry, level1, level2, level3] - not taken from anything, just imposed.
entrants = [0, 1, 2, 3]
start = 2;
neighbors_start = 108;
fields = 6;


for y in 1:size(yearins)[1]
	market_start = yearins[y][2]
	market_end = yearins[y][3]
#	market_frame = dataf[market_start:market_end, :]
	mkt_fips = yearins[y][1]
		for year in yearins[y][4:end]
			println(year)
			# There may be times when this is empty - check that.
	#		year_frame = dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year), :]
			level1 = dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:level1_hospitals0][1]
			level2 = dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:level2solo_hospitals0][1]
			level3 = dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:level3_hospitals0][1]
			fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid])) # this is run again below to catch entrants each iteration
			# What do I want to track over the whole history? fid, solo state, int state, action chosen, probability of choice, demand.  Aggregate: prob, levels.
			# Also think forward: demand realized.
			state_history = [zeros(1, fields*size(fids)[1]) 1 level1 level2 level3; zeros(T, fields*(size(fids)[1]) + 4)]
			# Includes values for the initial market configuration
			for n in 1:size(fids)[1]
				el = fids[n]
				a = ((dataf[:,:fid].==el)&(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year))
				state_history[1, (n-1)*fields + 1] = el # Change
				state_history[1, (n-1)*fields + 2] = dataf[a,:act_int][1]
				state_history[1, (n-1)*fields + 3] = dataf[a,:act_solo][1]
				state_history[1, (n-1)*fields + 4] = 1 #probability is 1 for the first action
				state_history[1, (n-1)*fields + 5] = 10 # No action at the first period?  Or should it be 10?
				#state_history[1, (n-1)*fields + 6] = # whatever demand is
			end
			# Record aggregte initial values
			state_history[1, (size(fids)[1])*fields+1] = level1 ;
			state_history[1, (size(fids)[1])*fields+2] = level2 ;
			state_history[1, (size(fids)[1])*fields+3] = level3 ;
			state_history[1, (size(fids)[1])*fields+4] = 1; # initial probability.
			# Start simulation here:
			for i = start:T+1
				if i%50 == 0
						println("Period ",i, " in Market-year ", year, " ", mkt_fips)
				end
				fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid])) # needs to be updated each round to catch entrants
						for fid in 1:size(fids)[1] # this has to be handled separately for each hospital, due to the geography issue
							el = fids[fid] # the dataframe is mutable.
							a = ((dataf[:,:fid].==el)&(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year))
							if sum(a) > 1
								println("two entries for ", el, " ", year, " ", mkt_fips)
							end
							if  ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (0,0)) # level 1, actions:
								probs1 = logitest((0,0), level1, level2, level3, convert(Array, [dataf[a,:lev105][1]; dataf[a,:lev205][1]; dataf[a,:lev305][1]; dataf[a,:lev1515][1]; dataf[a,:lev2515][1]; dataf[a,:lev3515][1]; dataf[a,:lev11525][1]; dataf[a,:lev21525][1]; dataf[a,:lev31525][1]]) )
								if probs1 == ValueException
									println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
									break
								end
								# all_hosp_probs[fid] = hcat(el, year_frame[a, :act_int], year_frame[a,:act_solo], 10, probs[1], 2, probs[2], 1, probs[3], 11, probs[4])
								# Reassign action choices -
								dataf[a, :choicenum0] = 10
								dataf[a, :pr_ch_0] = probs1[1]
								dataf[a, :choicenum1] = 2
								dataf[a, :pr_ch_1] = probs1[2]
								dataf[a, :choicenum2] = 1
								dataf[a, :pr_ch_2] = probs1[3]
								dataf[a, :choicenum3] = 11
								dataf[a, :pr_ch_3] = probs1[4]
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
									dataf[a,:act_int] = 1
									dataf[a,:act_solo] = 0
								elseif action1 == 1
									chprob = probs1[3]
									dataf[a,:act_int] = 0
									dataf[a,:act_solo] = 1
								elseif action1 == 11
									chprob = probs1[4]
									dataf[a,:act_int] = -999
									dataf[a,:act_solo] = -999
								else
									println("Fail")
									println("Bad Action Chosen by", el)
									break
								end
								# Set own distance counts to 0 for all categories
								(dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
								# This is the logical place to recompute the demand.
								# write out state values - in blocks:
								state_history[i, (fid-1)*fields + 1] = el # Change
								state_history[i, (fid-1)*fields + 2] = dataf[a,:act_int][1]
								state_history[i, (fid-1)*fields + 3] = dataf[a,:act_solo][1]
								state_history[i, (fid-1)*fields + 4] = chprob
								state_history[i, (fid-1)*fields + 5] = action1
								#state_history[i, (fid-1)*fields + 6] = demand
							elseif ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (1,0)) #level 2, actions:
								# Can't evaluation as nexti here if they are initialized to negative 1
								probs2 = logitest((1,0), level1, level2, level3, convert(Array, [dataf[a,:lev105][1]; dataf[a,:lev205][1]; dataf[a,:lev305][1]; dataf[a,:lev1515][1]; dataf[a,:lev2515][1]; dataf[a,:lev3515][1]; dataf[a,:lev11525][1]; dataf[a,:lev21525][1]; dataf[a,:lev31525][1]]) )
								if probs2 == ValueException
									println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
									break
								end
								dataf[a, :choicenum0] = 5
								dataf[a, :pr_ch_0] = probs2[1]
								dataf[a, :choicenum1] = 10
								dataf[a, :pr_ch_1] = probs2[2]
								dataf[a, :choicenum2] = 6
								dataf[a, :pr_ch_2] = probs2[3]
								dataf[a, :choicenum3] = 11
								dataf[a, :pr_ch_3] = probs2[4]
								# Action:
								action2 = sample([5, 10, 6, 11] ,WeightVec([probs2[1], probs2[2], probs2[3], probs2[4]]))
								if action2 == 5
									chprob = probs2[1]
									dataf[a,:act_int] = 0
									dataf[a,:act_solo] = 0
								elseif action2 == 10
									chprob = probs2[2]
									# no change to state
								elseif action2 == 6
									chprob = probs2[3]
									dataf[a,:act_int] = 1
									dataf[a,:act_solo] = 0
								elseif action2 == 11
									chprob = probs2[4]
									dataf[a,:act_int] = -999
									dataf[a,:act_solo] = -999
								else
									println("Fail")
									println("Bad Action Chosen by", el)
									break
								end
								# Set own distance counts to 0 for all categories
								(dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
								# write out state values - in blocks:
								state_history[i, (fid-1)*fields + 1] = el
								state_history[i, (fid-1)*fields + 2] = dataf[a,:act_int][1]
								state_history[i, (fid-1)*fields + 3] = dataf[a,:act_solo][1]
								state_history[i, (fid-1)*fields + 4] = chprob
								state_history[i, (fid-1)*fields + 5] = action2
								#state_history[i, (fid-1)*fields + 6] = demand
							elseif  ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (0,1)) #level 3, actions:
								probs3 = logitest((0,1), level1, level2, level3, convert(Array, [dataf[a,:lev105][1]; dataf[a,:lev205][1]; dataf[a,:lev305][1]; dataf[a,:lev1515][1]; dataf[a,:lev2515][1]; dataf[a,:lev3515][1]; dataf[a,:lev11525][1]; dataf[a,:lev21525][1]; dataf[a,:lev31525][1]]) )
								if probs3 == ValueException
									println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
									break
								end
								dataf[a, :choicenum0] = 4
								dataf[a, :pr_ch_0] = probs3[1]
								dataf[a, :choicenum1] = 3
								dataf[a, :pr_ch_1] = probs3[2]
								dataf[a, :choicenum2] = 10
								dataf[a, :pr_ch_2] = probs3[3]
								dataf[a, :choicenum3] = 11
								dataf[a, :pr_ch_3] = probs3[4]
								# Action:
								action3 = sample([4, 3, 10, 11] ,WeightVec([probs3[1], probs3[2], probs3[3], probs3[4]]))
								if action3 == 3
									chprob = probs3[2]
									dataf[a, :act_int] = 1
									dataf[a, :act_solo] = 0
								elseif action3 == 4
									chprob = probs3[1]
									dataf[a, :act_int] = 0
									dataf[a, :act_solo] =0
								elseif action3 == 10
									chprob = probs3[3]
									# no change to state
								elseif action3 == 11
									chprob = probs3[4]
									dataf[a, :act_int] = -999
									dataf[a, :act_solo] = -999
								else
									println("Fail")
									println("Bad Action Chosen by", el)
									break
								end
								# Set own distance counts to 0 for all categories
								(dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
								# write out state values - in blocks:
								state_history[i, (fid-1)*fields + 1] = el
								state_history[i, (fid-1)*fields + 2] = dataf[a,:act_int][1]
								state_history[i, (fid-1)*fields + 3] = dataf[a,:act_solo][1]
								state_history[i, (fid-1)*fields + 4] = chprob
								state_history[i, (fid-1)*fields + 5] = action3
								#state_history[i, (fid-1)*fields + 6] = demand
							elseif ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (-999,-999)) # has exited.
								# No new actions to compute, but record.
								state_history[i, (fid-1)*fields + 1] = el
								state_history[i, (fid-1)*fields + 2] = dataf[a,:act_int][1]
								state_history[i, (fid-1)*fields + 3] = dataf[a,:act_solo][1]
								state_history[i, (fid-1)*fields + 4] = 1 # exit is absorbing, so the choice prob is always 1
								state_history[i, (fid-1)*fields + 5] = 0 # no action is taken.
								# Set own distance counts to 0 for all categories
								(dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
							end
						end
						# Count facilities by distance and map results to neighboring hospitals
						# here the issue is that, for hospitals in neighboring counties, we haven't set the levX_YZ values to 0
						for f in fids
							own_fac = (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :act_solo][1], dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :act_int][1])
							for j = neighbors_start:(2):size(dataf)[2]
								if (!isna(dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j][1])) & (!isna(dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]))
									if (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), j+1][1] >0) & (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]< 5)
										if own_fac == (0,0)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev105] += 1
										elseif own_fac == (1,0)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev205] += 1
										elseif own_fac == (0,1)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev305] += 1
										elseif own_fac == (-999,-999)
											# do nothing - firm exited.
										else
											println("Bad Own Facility Code", f, j, own_fac)
										end
									elseif (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), j+1][1] >5) & (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]< 15)
										if own_fac == (0,0)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev1515] += 1
										elseif own_fac == (1,0)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev2515] += 1
										elseif own_fac == (0,1)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev3515] += 1
										elseif own_fac == (-999,-999)
											# do nothing - firm exited.
										else
											println("Bad Own Facility Code", f, j, own_fac)
										end
									elseif (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]> 15) & (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]< 25)
										if own_fac == (0,0)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev11525] += 1
										elseif own_fac == (1,0)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev21525] += 1
										elseif own_fac == (0,1)
											dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev31525] += 1
										elseif own_fac == (-999,-999)
											# do nothing - firm exited.
										else
											println("Bad Own Facility Code", f, j, own_fac)
										end
									else
										println("Bad Distance at frame ", f, j)
									end
								end
							end
						end
				# Sum the levels for next period:
				level1 = 0; level2 = 0; level3 = 0;
				update_mkt = ((dataf[:year].==year)&(dataf[:fipscode].==mkt_fips)&(dataf[:act_int].!=-999)&(dataf[:act_solo].!=-999))
				total = sum(update_mkt)
				intens = sum(dataf[update_mkt, :act_int])
				solo = sum(dataf[update_mkt, :act_solo])
				nones = sum(update_mkt) - intens - solo
						if (nones < 0) | (intens > total) | (solo > total) | (nones + intens + solo != total)
							println("Bad market size computations")
							println("Total ", total, " Level1 ", nones, " Level2 ", solo, " Level 3 ", intens, " Fips: ", mkt_fips, " Year ", year)
						else
								level1 = nones
								level2 = solo
								level3 = intens
						end

						# Entry draw
				entrypairs = hcat(entrants, entryprobs)
				entrantsp = WeightVec(entryprobs)
				newentrant = sample(entrants, entrantsp)
				entrantout = [newentrant, entrypairs[findfirst(entrypairs[:,1], newentrant), 2]]'
				b = (dataf[:fipscode].== mkt_fips)&(dataf[:year].==year) # Market-year observations, whether exited or not.
						if newentrant> 0
							# Eventually fix the fact that the neighbors here are not going to be exactly right.
							# I need to add the fact that the entrants are not being recorded in the "neighbors" section so no one is counting distances to them.
							ent_lat = mean(dropna(dataf[b,:v15])) + rand(Normal(0, 0.1), 1) # 0.1 degrees latitude should be about 6-7 miles.
							ent_lon = mean(dropna(dataf[b,:v16])) + rand(Normal(0, 0.1), 1)
							newrow = dataf[b,:][1,:] # create new dataframe row, duplicating existing.  Takes first row of current
							newrow[:facility] = convert(UTF8String, "Entrant $y $year")
							newrow[:fid] = sample(collect((maximum(dataf[:fid])+1):(maximum(dataf[:fid])+5))) # new facility will always have largest fid
							newrow[:id] = sample(collect((maximum(dataf[:id])+1):(maximum(dataf[:id])+5)))
							newrow[:fipscode] = mkt_fips
							newrow[:location] = convert(UTF8String, "entrant - see v15 v16")
							newrow[:city] = convert(UTF8String, "Entrant - unspecified")
							newrow[:firstyear] = year
							newrow[:v15] = ent_lat
							newrow[:v16] = ent_lon
							if newentrant == 1
								level1 += 1
								newrow[:act_int] = 0
								newrow[:act_solo] = 0
							elseif newentrant == 2
								level2 += 1
								newrow[:act_int] = 0
								newrow[:act_solo] = 1
							elseif newentrant == 3
								level3 += 1
								newrow[:act_int] = 0
								newrow[:act_solo] = 1
							end
							# Handle appending these entrants to the (neighbor, distance) section
							# need to check all of the other fids in the market-year (in b)
							for row in eachrow(dataf[b,:])
								if  (distance(ent_lat[1], ent_lon[1], row[:v15], row[:v16]) < 25)
									count = 0 # only want to make an entry once - this is a dumb way
									for c in neighbors_start:(2):size(dataf)[2]
										if (isna(row[c]))&(count == 0)
											# to make changes I need to search for the row in the original DF matching these characteristics.
											dataf[(dataf[:fid].==row[:fid])&(dataf[:id].==row[:id])&(dataf[:fipscode].==row[:fipscode]), c ]= newrow[:fid]
											dataf[(dataf[:fid].==row[:fid])&(dataf[:id].==row[:id])&(dataf[:fipscode].==row[:fipscode]), c+1]= distance(ent_lat[1], ent_lon[1], row[:v15], row[:v16])
											count += 1
										end
									end
								end
							end
							# Append the neighbors to the new entrant's frame too
							for elem in 1:size(fids)[1]
								neighb = fids[elem]
								neighb_lat = dataf[(dataf[:fid].==neighb)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :v15][1]
								neighb_lon = dataf[(dataf[:fid].==neighb)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :v16][1]
								td = distance(ent_lat[1], ent_lon[1], neighb_lat, neighb_lon)
								if  td < 25
									newrow[neighbors_start+2*(elem-1)] = neighb # fid
									newrow[neighbors_start+2*(elem-1)+1] = td
									if (td > 0) & (td < 5)
										if (newrow[:act_solo][1], newrow[:act_int][1]) == (0,0)
											newrow[:lev105] += 1
										elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (1,0)
											newrow[:lev205] += 1
										elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (0,1)
											newrow[:lev305] += 1
										else
											println("Bad facility in Entrant 1")
										end
									elseif (td > 5) & (td < 15)
										if (newrow[:act_solo][1], newrow[:act_int][1]) == (0,0)
											newrow[:lev1515] += 1
										elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (1,0)
											newrow[:lev2515] += 1
										elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (0,1)
											newrow[:lev3515] += 1
										else
											println("Bad facility in Entrant 2")
										end
									elseif (td > 15) & (td < 25)
										if (newrow[:act_solo][1], newrow[:act_int][1]) == (0,0)
											newrow[:lev11525] += 1
										elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (1,0)
											newrow[:lev21525] += 1
										elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (0,1)
											newrow[:lev31525] += 1
										else
											println("Bad facility in Entrant 3")
										end
									else
										println("Bad distance measured from entrant")
									end
								end
							end
							# Add the new record to the dataframe.
							append!(dataf, newrow)
							# append value to fids
							push!(fids, newrow[:fid][1])
							# Reshape state history: fid, solo state, int state, action chosen, probability of choice, demand. [newrow[:fid], 999, 999, 0, 1, 0]
							state_history = vcat(hcat(state_history[1:i,1:end-4], repmat([newrow[:fid][1] 999 999 0 1 0], i, 1), state_history[1:i, end-3:end]), zeros((T-i+1), size(fids)[1]*fields+4 ))
						end
				# Aggregate Probability of Action:
				tprob = 1
				for els in 4:fields:(size(state_history[i,:])[2]-4) # 4 is the relevant column
					if (state_history[i,els] < 0) | (state_history[i,els]>1)
 							println("Bad probability at row ", i, " ", els )
					else
						tprob = tprob*state_history[i,els]
					end
				end
				# Accounting for entry in aggregate prob of action -
				if newentrant >0
					if newentrant == 1
						tprob = tprob*entryprobs[2]
					elseif newentrant == 2
						tprob = tprob*entryprobs[3]
					elseif newentrant == 3
						tprob = tprob*entryprobs[4]
					end
				end
				state_history[i, (size(fids)[1])*fields+4] = state_history[i-1, (size(fids)[1])*fields+4]*tprob  #prob of ending up at previous state * current transition prob

				# Total number of firms
				state_history[i, (size(fids)[1])*fields+1] = level1 ;
				state_history[i, (size(fids)[1])*fields+2] = level2 ;
				state_history[i, (size(fids)[1])*fields+3] = level3 ;
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
