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
# using DataFrames
# using DataArrays
# using Distributions
#
# # Include necessary functions
# # include("/Users/austinbean/Desktop/dynhosp/combgen.jl")
# # include("/Users/austinbean/Desktop/dynhosp/nckr.jl")
# # include("/Users/austinbean/Desktop/dynhosp/probfinder.jl")
# # include("/Users/austinbean/Desktop/dynhosp/probfind2.jl")
# # include("/Users/austinbean/Desktop/dynhosp/tuplefinder.jl")
# include("/Users/austinbean/Desktop/dynhosp/LogitEst.jl")
# include("/Users/austinbean/Desktop/dynhosp/Distance.jl")
# include("/Users/austinbean/Desktop/dynhosp/Simulator.jl")
# include("/Users/austinbean/Desktop/dynhosp/PerturbSimulation.jl")
# include("/Users/austinbean/Desktop/dynhosp/DynamicValue.jl")
# include("/Users/austinbean/Desktop/dynhosp/DemandModel.jl")
#
#
# # Import Data
# data1 = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
# notmissing = findin(isna(data1[:fipscode]), false);
# data1 = data1[notmissing, :];
#
# regcoeffs = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Choice Model.csv", header = true);
#
# # Individual level demands -
# # DO NOT CHANGE THE NAME "people" - it will mess up fidfinder in DemandFunction.jl
# people = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 Individual Choices.csv", header = true);
#
# # Check the NFP status variable in the above
# # Coefficients on the demand model:
# modcoeffs = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 Model.csv", header = true);
# distance_c = modcoeffs[1, 2]
# distsq_c = modcoeffs[2, 2]
# neoint_c = modcoeffs[3, 2]
# soloint_c = modcoeffs[4, 2]
# closest_c = modcoeffs[5, 2]
# distbed_c = modcoeffs[6, 2]
#
# demandmodelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]
#

# For use in the demand model::
type  RowSizeError  <: Exception end



# Enumerate all of the types::
# a = Set()
# for el in people.columns
#   push!(a, typeof(el))
# end
#
# # This is needed to clean out the missing values among fids.  Changes them to 0.
#
# # Think about one change - if demand model coeff is positive, change value to large negative,
# # if negative, change to large positive.  Then it's basically impossible for that to be the choice.
#
#
for i in names(people)
  if ( typeof(people[i]) == DataArrays.DataArray{Float64,1} )
    people[DataFrames.isna(people[i]), i] = 0
  elseif (typeof(people[i]) == DataArrays.DataArray{Int64,1})
    people[DataFrames.isna(people[i]), i] = 0
  elseif typeof(people[i]) == DataArrays.DataArray{ByteString,1}
    # A dumb way to make sure no one chooses a missing facility: set covariate values to large numbers
    # with opposite signs of the corresponding coefficients from modelparameters.
    # This does that by looking at missing NAMES, not fids.
    people[DataFrames.isna(people[i]), people.colindex.lookup[i]+2] = -sign(neoint_c)*999
    people[DataFrames.isna(people[i]), people.colindex.lookup[i]+8] = -sign(soloint_c)*999
    people[DataFrames.isna(people[i]), i] = "NONE"
  end
  if sum(size(people[DataFrames.isna(people[i]), i]))>0
    print(i, "\n")
  end
end



### Collect Basic Information ###

# locates starting and ending points of all separate fipscode indices
allindices = [ (x, findfirst(data[:,fipscodeloc], x), findlast(data[:,fipscodeloc], x)) for x in unique(data[:,fipscodeloc]) ];
# Next one also works in case I decide an array of arrays is better than an array of tuples
#indices = [ [x, findfirst(df[:fipscode], x), findlast(df[:fipscode], x)] for x in unique(df[:fipscode]) ]
# Next one stores all the years appearing in each fips code
yearins = [ [x; findfirst(data[:,fipscodeloc], x); findlast(data[:,fipscodeloc], x ); unique( data[findfirst(data[:,fipscodeloc], x):findlast(data[:,fipscodeloc], x ) , yearloc]  ) ] for x in unique(data[:,fipscodeloc])  ]




# Parameters of the shock distribution
# This is wrong - I don't draw shocks from this, I draw conditional means,
# conditional on being the maximum.
# dist_μ = 0;
# dist_σ = 1;
# dist_ξ = 0;
# srand(123)
# d = GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# # I do need the constant:
# γ = eulergamma;


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

# β = 0.95;
# T = 100;
# choices = 4;
# α₂ = 0.07;  # Fraction of patients admitted to NICU lev 2 on average (PA Data)
# α₃ = 0.13; # Fraction of patients admitted to NICU Lev 3 on average (PA Data)
# # Actual entry probabilities will be substituted in later.
# #entryprobs = [0.9895, 0.008, 0.0005, 0.002] # [No entry, level1, level2, level3] - not taken from anything, just imposed.
# #entryprobs = [1.0, 0.0, 0.0, 0.0] # back-up entry probs with no entry for faster work.
# entrants = [0, 1, 2, 3]
# sim_start = 2;
# neighbors_start = 108;
# fields = 7;

# I should have 356 FIDs x 22 years of hospitals.

#
# mkt_fips = 48001;
# year = 2003;
# dataf = deepcopy(data1);
# fids = convert(Array, sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid])))
# peoples = people[fidfinder(fids, people, "people"),:];
# peoplesub = deepcopy(peoples);
# Mainfun(dataf, peoplesub, 48001, 2003, demandmodelparameters, entryprobs, fids)


#=
Future direction - write this to start where it left off.
# Start with some smaller markets first:

=#
monopoly = Array{Int64}(0)
duopoly = Array{Int64}(0)
triopoly = Array{Int64}(0)
tetrapoly = Array{Int64}(0)
nopoly = Array{Int64}(0)

for el in yearins
  unqfids = [x for x in unique(data1[el[2]:el[3],fidloc]).data]
  if size(unqfids)[1] == 1
    print("Fipscode Monopoly: ", unique(data1[el[2]:el[3], :fipscode]), "\n")
    push!(monopoly, data1[el[3], :fipscode])
  elseif size(unqfids)[1] == 2
    print("Fipscode Duopoly: ", unique(data1[el[2]:el[3], :fipscode]), "\n")
    push!(duopoly, data1[el[3], :fipscode])
  elseif size(unqfids)[1] == 3
    print("Fipscode Triopoly: ", unique(data1[el[2]:el[3], :fipscode]), "\n")
    push!(triopoly, data1[el[3], :fipscode])
  elseif size(unqfids)[1] == 4
    print("Fipscode Tetrapoly: ", unique(data1[el[2]:el[3], :fipscode]), "\n")
    push!(tetrapoly, data1[el[3], :fipscode])
  elseif size(unqfids)[1] > 4
    print("Fipscode N-opoly: ", unique(data1[el[2]:el[3], :fipscode]), "\n")
    print("Fipscode Hospitals: ", size(unqfids)[1], "\n")
    push!(nopoly, data1[el[3], :fipscode])
  end
end


container = zeros(1, 183)

# Open the existing saved data:
fout1 = readtable("/Users/austinbean/Desktop/dynhosp/simulationresults.csv")
donefips  = [x for x in unique(fout1[:fipscode])]


entryprobs = [0.9895, 0.008, 0.0005, 0.002]



for y in 1:size(yearins,1)    #size(duopoly)[1]
    mkt_fips = yearins[y][1] #duopoly[y]
    if !(in(mkt_fips, donefips)) # this is going to do new fipscodes only
      print("Market FIPS Code ", mkt_fips, "\n")
      	for year in [ 2003 2004 2005 2006]   #yearins[y][4:end] # can do all years or several.
          dataf = deepcopy(data);
          fids = convert(Array, sort!(unique(dataf[(dataf[:,fipscodeloc].==mkt_fips)&(dataf[:, yearloc].==year),fidloc])))
          numfids = size(fids)[1]
          peoples = people[fidfinder(fids, people),:];
          global peoplesub # the function below doesn't see "peoplesub" due to scope rules, unless it is declared as a global.
          peoplesub = deepcopy(peoples);
          # print("exists?: ", size(peoplesub), "\n")
          container = [container; Mainfun(dataf, peoplesub, mkt_fips, year, demandmodelparameters, entryprobs, fids)]
      end
    end
end


# Add column names to the new data:

colnames = Array{Symbol}(:0)
push!(colnames, :fipscode)
push!(colnames, :fid)
push!(colnames, :year)
for elem in ["EQ", "NEQ"]
  for j = 1:3
    for k = 0:25
      push!(colnames, parse("$elem"*"Lev$j"*"Comp$k"))
    end
  end
  for x in [1 2 3]
    for y in [1 2 3 "EX"]
      if x != y
        push!(colnames, parse("$elem"*"Trans$x$y"))
      end
    end
  end
  push!(colnames, parse("$elem"*"Enter1"))
  push!(colnames, parse("$elem"*"Enter2"))
  push!(colnames, parse("$elem"*"Enter3"))
end

output1 = convert(DataFrame, container);

names!(output1, colnames)
output1 = output1[ output1[:fipscode].>0 ,:]

# Append new data to the existing data:
append!(fout1, output1)
writetable("/Users/austinbean/Desktop/dynhosp/simulationresults.csv", output1)
















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
