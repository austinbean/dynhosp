# Module for Main Work:
# push this change to the LOAD_PATH: push!(LOAD_PATH, "/Users/austinbean/Desktop/dynhosp")
#
# using DataFrames
# @everywhere using DataFrames
#
# using Distributions
# @everywhere using Distributions

module ProjectModule

  using DataFrames
#  using DataArrays
  using Distributions

  include("/Users/austinbean/Desktop/dynhosp/LogitEst.jl")
  include("/Users/austinbean/Desktop/dynhosp/Distance.jl")
  include("/Users/austinbean/Desktop/dynhosp/Simulator.jl")
  include("/Users/austinbean/Desktop/dynhosp/PerturbSimulation.jl")
  include("/Users/austinbean/Desktop/dynhosp/DynamicValue.jl")
  include("/Users/austinbean/Desktop/dynhosp/DemandModel.jl")
  include("/Users/austinbean/Desktop/dynhosp/MainFunction.jl")
  include("/Users/austinbean/Desktop/dynhosp/PerturbAction.jl")
  include("/Users/austinbean/Desktop/dynhosp/ParallelFunctions.jl")


  export LogitEst, Simulator, PerturbSimulator, logitest, Mainfun, DemandModel, EntrantsU, rowchange, rowfindfid, fidfinder, sendto, passobj


end # of module


#=
neighbors_start = 108;
entryprobs = [0.99, 0.004, 0.001, 0.005] # [No entry, level1, level2, level3] - not taken from anything, just imposed.
entrants = [0, 1, 2, 3]
fields = 7; # if fields updated, update reshaping of state history
T = 100;

yearins = [ [x; findfirst(dataf[:fipscode], x); findlast(dataf[:fipscode], x ); unique( dataf[findfirst(dataf[:fipscode], x):findlast(dataf[:fipscode], x ) , :year]  ) ] for x in unique(dataf[:fipscode])  ]
mkt_fips = yearins[10][1]
year = 2011
fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid]))



dataf = DataFrames.readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
notmissing = findin(DataFrames.isna(dataf[:fipscode]), false);
dataf = dataf[notmissing, :];

  # Handle missing dataframe elements::
  for i in names(dataf)
      if ( typeof(dataf[i]) == DataArrays.DataArray{Float64,1} )
          dataf[DataFrames.isna(dataf[i]), i] = 0
      elseif (typeof(dataf[i]) == DataArrays.DataArray{Int64,1})
          dataf[DataFrames.isna(dataf[i]), i] = 0
      elseif typeof(dataf[i]) == DataArrays.DataArray{ByteString,1}
          dataf[DataFrames.isna(dataf[i]), i] = "NONE"
    end
      if sum(size(dataf[DataFrames.isna(dataf[i]), i]))>0
      print(i, "\n")
    end
  end

  data = convert(Matrix, dataf);



# TESTING --- To run a test, replace the value in the first bracket in mkt_fips = yearins[ ][1], and choose a year.

# Import Data
dataf = DataFrames.readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
notmissing = findin(isna(dataf[:fipscode]), false);
dataf = dataf[notmissing, :];
yearins = [ [x; findfirst(dataf[:fipscode], x); findlast(dataf[:fipscode], x ); unique( dataf[findfirst(dataf[:fipscode], x):findlast(dataf[:fipscode], x ) , :year]  ) ] for x in unique(dataf[:fipscode])  ]
mkt_fips = yearins[10][1]
year = 2011
fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid]))

# Handle missing dataframe elements::
for i in names(dataf)
    if ( typeof(dataf[i]) == DataArrays.DataArray{Float64,1} )
        dataf[isna(dataf[i]), i] = 0
    elseif (typeof(dataf[i]) == DataArrays.DataArray{Int64,1})
        dataf[isna(dataf[i]), i] = 0
    elseif typeof(dataf[i]) == DataArrays.DataArray{ByteString,1}
        dataf[isna(dataf[i]), i] = "NONE"
  end
    if sum(size(dataf[isna(dataf[i]), i]))>0
    print(i, "\n")
  end
end

data = convert(Matrix, dataf)

# Now send to all workers:
sendto(lis; data=data)


# Load the people
people = DataFrames.readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);

# Enumerate all of the types::
a = Set()
for el in people.columns
  push!(a, typeof(el))
end

# This is needed to clean out the missing values among fids.  Changes them to 0.
for i in names(people)
  if typeof(people[i]) != DataArrays.DataArray{ByteString,1}
    people[DataFrames.isna(people[i]), i] = 0
  elseif typeof(people[i]) == DataArrays.DataArray{ByteString,1}
    people[DataFrames.isna(people[i]), i] = "NONE"
  end
end

peoples = convert(Matrix, people);


# Send to all workers:
sendto(lis; people=people)



modcoeffs = DataFrames.readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Model.csv", header = true);
distance_c = modcoeffs[1, 2]
distsq_c = modcoeffs[2, 2]
neoint_c = modcoeffs[3, 2]
soloint_c = modcoeffs[4, 2]
closest_c = modcoeffs[5, 2]
distbed_c = modcoeffs[6, 2]

demandmodelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]

sendto(lis, demandmodelparameters = demandmodelparameters)



fipscodeloc = dataf.colindex.lookup[:fipscode]
sendto(lis,  fipscodeloc=fipscodeloc )

yearloc = dataf.colindex.lookup[:year]
sendto(lis,  yearloc=yearloc )

level1_hospitals0loc = dataf.colindex.lookup[:level1_hospitals0]
sendto(lis,  level1_hospitals0loc=level1_hospitals0loc )

level2solo_hospitals0loc = dataf.colindex.lookup[:level2solo_hospitals0]
sendto(lis,  level2solo_hospitals0loc=level2solo_hospitals0loc )

level3_hospitals0loc = dataf.colindex.lookup[:level3_hospitals0]
sendto(lis,  level3_hospitals0loc=level3_hospitals0loc )

fidloc = dataf.colindex.lookup[:fid]
sendto(lis,  fidloc=fidloc )

act_intloc = dataf.colindex.lookup[:act_int]
sendto(lis,  act_intloc=act_intloc )

act_sololoc = dataf.colindex.lookup[:act_solo]
sendto(lis,  act_sololoc=act_sololoc )

lev105loc = dataf.colindex.lookup[:lev105]
sendto(lis,  lev105loc=lev105loc )

lev205loc = dataf.colindex.lookup[:lev205]
sendto(lis,  lev205loc=lev205loc )

lev305loc = dataf.colindex.lookup[:lev305]
sendto(lis,  lev305loc=lev305loc )

lev1515loc = dataf.colindex.lookup[:lev1515]
sendto(lis,  lev1515loc=lev1515loc )

lev2515loc = dataf.colindex.lookup[:lev2515]
sendto(lis,  lev2515loc=lev2515loc )

lev3515loc = dataf.colindex.lookup[:lev3515]
sendto(lis,  lev3515loc=lev3515loc )

lev11525loc = dataf.colindex.lookup[:lev11525]
sendto(lis,  lev11525loc=lev11525loc )

lev21525loc = dataf.colindex.lookup[:lev21525]
sendto(lis,  lev21525loc=lev21525loc )

lev31525loc = dataf.colindex.lookup[:lev31525]
sendto(lis,  lev31525loc=lev31525loc )

v15loc = dataf.colindex.lookup[:v15]
sendto(lis,  v15loc=v15loc )

v16loc = dataf.colindex.lookup[:v16]
sendto(lis,  v16loc=v16loc )

facilityloc = dataf.colindex.lookup[:facility]
sendto(lis,  facilityloc=facilityloc )

idloc = dataf.colindex.lookup[:id]
sendto(lis,  idloc=idloc )

locationloc = dataf.colindex.lookup[:location]
sendto(lis,  locationloc=locationloc )

cityloc = dataf.colindex.lookup[:city]
sendto(lis,  cityloc=cityloc )

firstyearloc = dataf.colindex.lookup[:firstyear]
sendto(lis,  firstyearloc=firstyearloc )


TotalBeds1loc = people.colindex.lookup[:TotalBeds1]
sendto(lis,  TotalBeds1loc=TotalBeds1loc )

TotalBeds2loc = people.colindex.lookup[:TotalBeds2]
sendto(lis,  TotalBeds2loc=TotalBeds2loc )


colparams = [fipscodeloc yearloc level1_hospitals0loc level2solo_hospitals0loc level3_hospitals0loc fidloc act_intloc act_sololoc lev105loc lev205loc lev305loc lev1515loc lev2515loc lev3515loc lev11525loc lev21525loc lev31525loc v15loc v16loc facilityloc idloc locationloc cityloc firstyearloc TotalBeds1loc TotalBeds2loc]


dist_μ = 0;
dist_σ = 1;
dist_ξ = 0;
srand(123)
d = GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# I do need the constant:
γ = eulergamma;


=#
