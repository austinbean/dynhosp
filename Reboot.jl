# reboot:
# include("/Users/austinbean/Desktop/dynhosp/Reboot.jl")


push!(LOAD_PATH, "/Users/austinbean/Desktop/dynhosp")

lis = addprocs()

using ProjectModule
using DataFrames
using Distributions


@everywhere begin
  # Import the hospital data and convert to a matrix -
    dataf = DataFrames.readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
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

  # Import the people and convert that data to a matrix
    people = DataFrames.readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);

    for i in names(people)
      if typeof(people[i]) != DataArrays.DataArray{ByteString,1}
        people[DataFrames.isna(people[i]), i] = 0
      elseif typeof(people[i]) == DataArrays.DataArray{ByteString,1}
        people[DataFrames.isna(people[i]), i] = "NONE"
      end
    end
    peoples = convert(Matrix, people);

  # Import the model coefficients
  modcoeffs = DataFrames.readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Model.csv", header = true);

  global const distance_c = modcoeffs[1, 2]
  global const distsq_c = modcoeffs[2, 2]
  global const neoint_c = modcoeffs[3, 2]
  global const soloint_c = modcoeffs[4, 2]
  global const closest_c = modcoeffs[5, 2]
  global const distbed_c = modcoeffs[6, 2]

  demandmodelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]

  global const idloc = dataf.colindex.lookup[:id]


end # of "begin" block



# Test functions:
a1 = Simulator(data, peoples, 2003, 48027, demandmodelparameters; T = 5) #tests logitest, DemandModel, distance
b1 = PerturbSimulator(data, peoples, 2003, 48027, demandmodelparameters, 273410; T = 5)
logitest( (1, 0), 2, 3, 1, zeros(9))
DemandModel(peoples, demandmodelparameters, [99999 1 0 120 32.96  -96.8385]) #implicitly tests EntrantsU as well
DynamicValue(b1, b1[1,1])

# Test parallel

p1 = remotecall_fetch(lis[1], DemandModel, peoples, demandmodelparameters, Array{Float64,2}())



# This is necessary for Simulator and PerturbSimulator - column numbers of data elements in "data" and "dataf"
#
# global const fipscodeloc = dataf.colindex.lookup[:fipscode]
# global const yearloc = dataf.colindex.lookup[:year]
# global const level1_hospitals0loc = dataf.colindex.lookup[:level1_hospitals0]
# global const level2solo_hospitals0loc = dataf.colindex.lookup[:level2solo_hospitals0]
# global const level3_hospitals0loc = dataf.colindex.lookup[:level3_hospitals0]
# global const fidloc = dataf.colindex.lookup[:fid]
# global const act_intloc = dataf.colindex.lookup[:act_int]
# global const act_sololoc = dataf.colindex.lookup[:act_solo]
# global const lev105loc = dataf.colindex.lookup[:lev105]
# global const lev205loc = dataf.colindex.lookup[:lev205]
# global const lev305loc = dataf.colindex.lookup[:lev305]
# global const lev1515loc = dataf.colindex.lookup[:lev1515]
# global const lev2515loc = dataf.colindex.lookup[:lev2515]
# global const lev3515loc = dataf.colindex.lookup[:lev3515]
# global const lev11525loc = dataf.colindex.lookup[:lev11525]
# global const lev21525loc = dataf.colindex.lookup[:lev21525]
# global const lev31525loc = dataf.colindex.lookup[:lev31525]
# global const v15loc = dataf.colindex.lookup[:v15]
# global const v16loc = dataf.colindex.lookup[:v16]
# global const facilityloc = dataf.colindex.lookup[:facility]

# global const locationloc = dataf.colindex.lookup[:location]
# global const cityloc = dataf.colindex.lookup[:city]
# global const firstyearloc = dataf.colindex.lookup[:firstyear]
#
# global const TotalBeds1loc = people.colindex.lookup[:TotalBeds1]
# global const TotalBeds2loc = people.colindex.lookup[:TotalBeds2]