
# Store the hospitals in data structures to reduce the complexity of writing out all of the values
# to matrices and keeping track of all of the indices, etc.

using Distributions
using DataFrames # WeightVec is in DataFrames, not Distributions.


push!(LOAD_PATH, "/Users/austinbean/Desktop/dynhosp")
push!(LOAD_PATH, "/dynhosp/dynhosp")
push!(LOAD_PATH, "/home/ubuntu/dynhosp/")
push!(LOAD_PATH, "/home1/04179/abean/dynhosp")
dir = pwd()
global pathdata = "";global pathpeople = "";global  pathprograms = "";
if dir == "/Users/austinbean/Desktop/dynhosp"
  global pathdata = "/Users/austinbean/Google Drive/Annual Surveys of Hospitals/"
  global pathpeople = "/Users/austinbean/Google Drive/Texas Inpatient Discharge/"
  global pathprograms = "/Users/austinbean/Desktop/dynhosp/"
elseif dir == "/dynhosp/dynhosp"
  global pathdata = dir*"/"
  global pathpeople = dir*"/"
  global pathprograms = dir*"/"
elseif (dir == "/home/ubuntu/Notebooks") | (dir == "/home/ubuntu/dynhosp")
  global pathdata = "/home/ubuntu/dynhosp/"
  global pathpeople = "/home/ubuntu/dynhosp/"
  global pathprograms = "/home/ubuntu/dynhosp/"
elseif (dir == "/home1/04179/abean/dynhosp")
  global pathdata = "/home1/04179/abean/dynhosp"
  global pathpeople = "/home1/04179/abean/dynhosp"
  global pathprograms = "/home1/04179/abean/dynhosp"
end
# Import the hospital data and convert to a matrix -
println("Importing Hosp Data")
  dataf = DataFrames.readtable(pathdata*"TX Transition Probabilities.csv", header = true);
  for i in names(dataf)
      if ( typeof(dataf[i]) == DataArrays.DataArray{Float64,1} )
          dataf[DataFrames.isna(dataf[i]), i] = 0
      elseif (typeof(dataf[i]) == DataArrays.DataArray{Int64,1})
          dataf[DataFrames.isna(dataf[i]), i] = 0
      elseif typeof(dataf[i]) == DataArrays.DataArray{ByteString,1}
          dataf[DataFrames.isna(dataf[i]), i] = "NONE"
      elseif typeof(dataf[i]) == DataArrays.DataArray{UTF8String,1}
            dataf[DataFrames.isna(dataf[i]), i] = "NONE"
    end
      if sum(size(dataf[DataFrames.isna(dataf[i]), i]))>0
      print(i, "\n")
    end
  end
  data = convert(Matrix, dataf);
#  dataf = 0; #set to zero to clear out.


type WTP
  w385::Vector
  w386::Vector
  w387::Vector
  w388::Vector
  w389::Vector
  w390::Vector
  w391::Vector
end

type DemandHistory
  demand385::Vector
  demand386::Vector
  demand387::Vector
  demand388::Vector
  demand389::Vector
  demand390::Vector
  demand391::Vector
end

type neighbors
  level105::Int
  level205::Int
  level305::Int
  level1515::Int
  level2515::Int
  level3515::Int
  level11525::Int
  level21525::Int
  level31525::Int
end

type hospital
  fid::Int64
  lat::Float64
  long::Float64
  name::AbstractString
  fipscode::Int
  level::Int
  levelhistory::Vector{Int}
  demandhist::DemandHistory
  wtphist::WTP
  chprobability::WeightVec
  # The logitest function takes the following:
  # logitest((0,0), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
  neigh::neighbors
  hood::Array{hospital, 1}
end


type Market
	config::Array{hospital, 1}
  collection::Dict{Int64, hospital} # create the dict with a comprehension to initialize
  fipscode::Int
end

type EntireState
  ms::Array{Market, 1}
  mkts::Dict{Int64, Market}   # Link markets by FIPS code via dictionary.
end

# Initialize Empty collection of markets:
Texas = EntireState(Array{hospital,1}(), Dict{Int64, Market}())
# See below for dictionary comprehension.

fips = unique(data[:,78])
# This adds a list of the markets covered to the whole state iterable.
for el in fips
  if el != 0
    el = eval(parse("m$el = Market( Array{hospital,1}(), Dict{Int64, hospital}() ,$el)"))
    push!(Texas.ms, el)
  end
end
# Write the market dictionary out as a comprehension:
Texas.mkts = [ m.fipscode => m for m in Texas.ms]

# fid - col 74, lat - col 94, long - col 95, name - col 82, fipscode - col 78, act_int - col 79, act_solo - col 80

data05 = data[(data[:,75].==2005), :] ;
lev105loc = 97; lev205loc = 98; lev305loc = 99; lev1515loc = 101; lev2515loc = 102; lev3515loc = 103; lev11525loc = 105; lev21525loc = 106; lev31525loc = 107;
for i = 1:size(data05,1)
  fips = data05[i, 78]
  if fips != 0
    level = 0
    if (data05[i, 79] == 1)&(data05[i,80]==0)
      level = 3
    elseif (data05[i, 79] == 0)&(data05[i,80]==1)
      level = 2
    else
      level = 1
    end
    push!(Texas.mkts[fips].config,
    hospital( data05[i, 74],
              data05[i,94],
              data05[i, 95],
              data05[i, 82],
              fips,
              level,
              [level],
              DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WeightVec([data[i,19], data[i,37], data[i,55], data[i, 73]]),
              neighbors(data[i, lev105loc], data[i,lev205loc ], data[i,lev305loc ], data[i,lev1515loc ], data[i,lev2515loc ], data[i, lev3515loc], data[i,lev11525loc ], data[i,lev21525loc ], data[i,lev31525loc]  ),
              Array{hospital,1}() ) )
  end
end


function MarketPrint(mkt::Market)
  for el in mkt.config
    println(el.name)
  end
end

function NeighborsPrint(mkt::Market)
  for el in mkt.config
    println(el.name, " ", el.neigh)
  end
end

function MktSize(n::neighbors, variety::Int)
  if variety == 1
    return n.lev105 + n.lev1515 + n.lev11525
  elseif variety == 2
    return n.lev205 + n.lev2515 + n.lev21525
  elseif variety == 3
    return n.lev305 + n.lev3515 + n.lev31525
  else
    return n.lev105 + n.lev1515 + n.lev11525 + n.lev205 + n.lev2515 + n.lev21525 + n.lev305 + n.lev3515 + n.lev31525
  end 
end

#=

# Test Hospitals:



a1 = hospital(123,
              31.20323,
              41.239453,
              "Hosp1",
              48001,
              1,
              Array{Int,1}(),
              DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WeightVec([0.1, 0.2]),
              Array{hospital,1}())


a2 = hospital(456,
              31.20323,
              41.239453,
              "Hosp2",
              48001,
              1,
              Array{Int,1}(),
              DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WeightVec([0.1, 0.2]),
              Array{hospital,1}())

a3 = hospital(789,
              31.20323,
              41.239453,
              "Hosp3",
              48001,
              1,
              Array{Int,1}(),
              DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WeightVec([0.1, 0.2]),
              Array{hospital,1}())


a4 = hospital(1011,
              31.20323,
              41.239453,
              "Hosp4",
              48001,
              1,
              Array{Int,1}(),
              DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
              WeightVec([0.1, 0.2]),
              Array{hospital,1}())


# Test Markets


m1 = Market([a1, a2], Dict{Int64, hospital}(), 48001)
m2 = Market([a3, a4], Dict{Int64, hospital}(), 48002)


# Initialize the above w/ empty dict, then create with:
# This should make it much, much easier to look everything up.
m1.collection = [ a.fid => a for a in m1.config]
m2.collection = [ a.fid => a for a in m2.config]

# Demonstrations:


# This works, since config::Array{hospital}
for el in m1.config
  println(el)
  println(el.demandhist.demand385)
end
# to add demand of some kind or other
for el in m1.config
  push!(el.demandhist.demand385, 210)
end


# Illustrations of iteration -

# Note that here the iterable object inside the type must be accessed to iterate.
for market in Texas.ms
  for hosp in market.config
    println(hosp.name)
  end
end


=#

#=
# These are *not* necessary.
# Because config is a kind of array, we can already iterate over the elements.
Base.start(M::Market) = M.config[1]
Base.next(M::Market, state) = (M.config[state], state+1)
Base.done(M::Market, state) = state > size(M.config, 1)
=#
