# Main

#=
Version Information:
Julia Version 0.4.3
Commit a2f713d (2016-01-12 21:37 UTC)
Platform Info:
  System: Darwin (x86_64-apple-darwin13.4.0)
  CPU: Intel(R) Core(TM) i7-5557U CPU @ 3.10GHz
  WORD_SIZE: 64
  BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
  LAPACK: libopenblas64_
  LIBM: libopenlibm
  LLVM: libLLVM-3.3

=#

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
include("/Users/austinbean/Desktop/dynhosp/DynamicValue.jl")
include("/Users/austinbean/Desktop/dynhosp/DemandModel.jl")


# Import Data
data1 = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
notmissing = findin(isna(data1[:fipscode]), false);
dataf = data1[notmissing, :];

regcoeffs = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Choice Model.csv", header = true);

# Individual level demands -
# DO NOT CHANGE THE NAME "people" - it will mess up fidfinder in DemandFunction.jl
people = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);

# Check the NFP status variable in the above
# Coefficients on the demand model:
modcoeffs = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Model.csv", header = true);
distance_c = modcoeffs[1, 2]
distsq_c = modcoeffs[2, 2]
neoint_c = modcoeffs[3, 2]
soloint_c = modcoeffs[4, 2]
closest_c = modcoeffs[5, 2]
distbed_c = modcoeffs[6, 2]

modelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]


# For use in the demand model::
type  RowSizeError  <: Exception end



# Enumerate all of the types::
a = Set()
for el in people.columns
  push!(a, typeof(el))
end

# This is needed to clean out the missing values among fids.  Changes them to 0.
for i in names(people)
  if typeof(people[i]) != DataArrays.DataArray{UTF8String,1}
    people[isna(people[i]), i] = 0
  end
end




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


- Must turn these vectors into valuations.

- Demand system.

=#

container = zeros(5000, 182)

for y in 1:size(yearins)[1]
	mkt_fips = yearins[y][1]
  print("Market FIPS Code ", mkt_fips)
		for year in yearins[y][4:end]
      print("Year ", year)
			fids = convert(Array, sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid]))) # returns a dataframe unless converted
			# Find the subset of people with those fids::
			# DO NOT CHANGE THIS NAME! DemandModel function will be screwed up!
			peoples = people[fidfinder(fids, people, "people"),:] # DO NOT CHANGE NAME
      # The values are going to get mangled, so need to copy them.
      peoplesub = deepcopy(peoples)
      dataf = deepcopy(data1)

				#Arguments: Simulator(dataf::DataFrame, peoplesub::DataFrame, year::Int64, mkt_fips::Int64,  state_history::Array{Float64,2}, demandmodelparameters::Array{Float64, 2}; T = 100, start = 2)
        print("Equilibrium Simulation, ", mkt_fips, " ", year, " ")
      states = Simulator(dataf, peoplesub, "peoplesub", year, mkt_fips, modelparameters, T = 100, sim_start = 2)

			# Non-equilibrium Play -
			# Entrants in dataframe now tagged with negative ID's.  Remake to remove them:
			dataf = dataf[dataf[:id].>= 0, :]

			for f in 1:size(fids)[1]
        pfid = fids[f]
  			#Arguments: function PerturbSimulator(dataf::DataFrame, peoplesub::DataFrame, subname::ASCIIString, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2}, pfid::Int64; disturb = 0.05, T = 100, sim_start = 2)
        print("Non - Equilibrium Simulation, ", mkt_fips, " ", year, " ")
        perturbed_history = PerturbSimulator(dataf, peoplesub, "peoplesub", year, mkt_fips, modelparameters, pfid, disturb = 0.01, T = 100, sim_start = 2)

  			# Here apply DynamicValue to the result of the simulations
  			# DynamicValue(state_history::Array, fac_fid::Float64; pat_types = 1, β = 0.95, T = 100, max_hosp = 25)
  			# output is in format: facility changes record, per-period visits record.
  			pfid_f = convert(Float64, pfid)
  			eq_change, eq_val  = DynamicValue(states, pfid_f; pat_types = 1, β = 0.95, T = 100, max_hosp = 25)
  			neq_change, neq_val = DynamicValue(perturbed_history, pfid_f; pat_types = 1, β = 0.95, T = 100, max_hosp = 25)
        i = y*(year-1989)
        container[i,:] = [pfid_f year eq_val eq_change neq_val neq_change]
        # Abandon Entrants again.
        dataf = dataf[dataf[:id].>= 0, :]
      end
		end
end


# Extract the count of visits to various states in this section
# Lev 1: (0,0)
# Lev 2: (1,0)
# Lev 3: (0,1)









###



#=

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


=#
