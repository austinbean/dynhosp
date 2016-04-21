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
include("/Users/austinbean/Desktop/dynhosp/Simulator.jl")
include("/Users/austinbean/Desktop/dynhosp/PerturbSimulation.jl")



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
• NB: 48029 is Bexar County

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
sim_start = 2;
neighbors_start = 108;
fields = 7;

#=
To work on::
- The dataframe will be modified and so the number of hospitals might change
within the loop.  Either copy it or reload it between running the first and second.

- Must turn these vectors into valuations.

- Demand system.

=#


for y in 1:size(yearins)[1]
	mkt_fips = yearins[y][1]
		for year in yearins[y][4:end]
			fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid]))

			# Equilibrium Play -
			state_history = [zeros(1, fields*size(fids)[1]) 1 0 0 0; zeros(T, fields*(size(fids)[1]) + 4)]
				#Arguments: Simulator(dataf::DataFrame, year::Int64, mkt_fips::Int64,  state_history::Array{Float64,2}; T = 100, start = 2)
			states = Simulator(dataf, year, mkt_fips, state_history, T = 100, sim_start = 2)

			# Non-equilibrium Play -
			# BE CAREFUL - DATAFRAME HAS CHANGED (why is that word purple?)
			# Entrants in dataframe now tagged with negative ID's.  Remake:
			dataf = dataf[dataf[:id].>= 0, :]

			pfid = fids[1]
			p_history = [zeros(1, fields*size(fids)[1]) 1 0 0 0; zeros(T, fields*(size(fids)[1]) + 4)]
				#Arguments: PerturbSimulator(dataf::DataFrame, year::Int64, mkt_fips::Int64,  state_history::Array{Float64,2}, pfid::Int64; disturb = 0.05, T = 100, sim_start = 2)
			perturbed_history = PerturbSimulator(dataf, year, mkt_fips, p_history, pfid, disturb = 0.01, T = 100, sim_start = 2)
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

	previous = state_history[row, 1] # records the current state as "previous for the next round"
end












###
