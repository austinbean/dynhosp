
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
  hood::Array{Int64, 1} # keep an array of fids here, rather than an array of hospitals, since the ref is circular.
end


type Market
	config::Array{hospital, 1}
  collection::Dict{Int64, hospital} # create the dict with a comprehension to initialize
  fipscode::Int
end

type EntireState
  ms::Array{Market, 1}
  mkts::Dict{Int64, Market}   # Link markets by FIPS code via dictionary.
  fipsdirectory::Dict{Int64,Int64} # Directory should be hospital fid / market fips
end

# Initialize Empty collection of markets:
Texas = EntireState(Array{hospital,1}(), Dict{Int64, Market}(), Dict{Int64, hospital}())
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
              Array{Int64,1}() ) )
  end
  # push all hospital fid/ fips pairs into the directory.
  Texas.fipsdirectory[data05[i, 74]] = fips # now for the whole state I can immediately figure out which market a hospital is in.
end

# Expand the market dictionaries so that they are filled with the hospitals
for el in Texas.ms
  el.collection = [ i.fid => i for i in el.config ]
  # I would like to append to each hospital a list of the others in the market.
  for hosp in el.config
    for hosp2 in el.config
      if hosp.fid != hosp2.fid
        push!(hosp.hood, hosp2.fid)
      end
    end
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

function NewEntrantLocation(mkt::Market)
# Takes the market, takes the mean location of all hospitals, adds normal noise to it.
  meanlat = 0
  meanlong = 0
  for el in mkt.config # over hospitals
    meanlat += el.lat
    meanlong += el.long
  end
  return [meanlat/size(mkt.config, 1) + rand(Normal(0, 0.1), 1)[1], meanlong/size(mkt.config, 1) + rand(Normal(0, 0.1), 1)[1]]
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

# When the facility level changes, the choices need to change too.
function ChoicesAvailable(h::hospital)
  if h.level == 1
    return [10 2 1 11]
  elseif h.level == 2
    return [5 10 6 11]
  elseif h.level == 3
    return [4 3 10 11]
  else # exited
    return [-999 -999 -999 -999]
  end
end

function LevelFunction(h::hospital, choice::Int64)
  # Takes a hospital record and a choice and returns the corresponding level.
  if h.level == 1
    if choice == 10
      return 1
    elseif choice == 2
      return 3
    elseif choice == 1
      return 2
    else # choice must be 11
      return -999
    end
  elseif h.level == 2
    if choice == 5
      return 1
    elseif choice == 10
      return 2
    elseif choice == 6
      return 3
    else # choice must be 11
      return -999
    end
  elseif h.level == 3
    if choice == 4
      return 1
    elseif choice == 3
      return 2
    elseif choice == 10
      return 3
    else # choice must be 11
      return -999
    end
  else # value must be -999
    return -999
  end
end


function NeighborAppend(elm::hospital, entrant::hospital)
  #=
  takes two hospital records, computes the distance between them and adds a 1 to the relevant record in the neighborhood type.
   And appends it to the hood of elm, which is a list of fids.
  =#
  dist = distance(elm.lat, elm.long, entrant.lat, entrant.long )
  if !in(entrant.fid, elm.hood)
    if dist < 25
      push!(elm.hood, entrant.fid)
      if dist<5
        if entrant.level == 1
          elm.neigh.level105 += 1
        elseif entrant.level == 2
          elm.neigh.level205 += 1
        elseif entrant.level == 3
          elm.neigh.level305 += 1
        end
      elseif (dist>5)&(dist<15)
        if entrant.level == 1
          elm.neigh.level1515 += 1
        elseif entrant.level == 2
          elm.neigh.level2515 += 1
        elseif entrant.level == 3
          elm.neigh.level3515 += 1
        end
      elseif (dist>15)
        if entrant.level == 1
          elm.neigh.level11525 += 1
        elseif entrant.level == 2
          elm.neigh.level21525 += 1
        elseif entrant.level == 3
          elm.neigh.level31525 += 1
        end
      end
    end
  end
end


function NeighborRemove(elm::hospital, entrant::hospital)
  #=
  takes two hospital records, computes the distance between them and adds a 1 to the relevant record in the neighborhood type.
   And appends it to the hood of elm, which is a list of fids.
  =#
  dist = distance(elm.lat, elm.long, entrant.lat, entrant.long )
  if in(entrant.fid, elm.hood)
    if dist < 25
      deleteat!(elm.hood, findin(elm.hood, entrant.fid))
      if dist<5
        if entrant.level == 1
          elm.neigh.level105 = max(elm.neigh.level105 -1, 0)
        elseif entrant.level == 2
          elm.neigh.level205 = max(elm.neigh.level205 -1, 0)
        elseif entrant.level == 3
          elm.neigh.level305 = max(elm.neigh.level305 -1, 0)
        end
      elseif (dist>5)&(dist<15)
        if entrant.level == 1
          elm.neigh.level1515 = max(elm.neigh.level1515 -1,0)
        elseif entrant.level == 2
          elm.neigh.level2515 = max(elm.neigh.level2515 -1,0)
        elseif entrant.level == 3
          elm.neigh.level3515 = max(elm.neigh.level3515 -1,0)
        end
      elseif (dist>15)
        if entrant.level == 1
          elm.neigh.level11525= max(elm.neigh.level11525 - 1,0)
        elseif entrant.level == 2
          elm.neigh.level21525= max(elm.neigh.level21525 - 1,0)
        elseif entrant.level == 3
          elm.neigh.level31525= max(elm.neigh.level31525 - 1,0)
        end
      end
    end
  end
end

function HospFindFirst(mkt::Market, hosp::hospital)
  for el in 1:size(mkt.config,1)
    if mkt.config[el].fid == hosp.fid
      return el
    end
  end
end

function MarketCleaner(mkt::Market)
  # Should take as an argument a whole market and then remove the records of any entrants.
  for el in mkt.config
    if el.fid < 0
      # Remove the hospital from the market array
      deleteat!(mkt.config, HospFindFirst(mkt, el)) # NB: HospFindFirst takes *market* as argument, but deleteat! takes *array*, i.e, market.config
      # Remove the hospital from the market dictionary
      pop!(mkt.collection, el.fid)
      # Remove from the list of neighbors in
      for others in mkt.config
        NeighborRemove(el, others)
      end
    end
  end
end


entrants = [0, 1, 2, 3]
entryprobs = [0.9895, 0.008, 0.0005, 0.002]
for i = 1:5
  for el in Texas.ms
    entrant = sample(entrants, WeightVec(entryprobs))
    if entrant != 0
      println("Entry! ", entrant , " FIPS ", el.fipscode)
      entloc = NewEntrantLocation(el) # called on the market
      newfid = -floor(rand()*1e6)
      entr = hospital( newfid,
                       entloc[1],
                       entloc[2],
                       " Entrant $newfid ",
                       el.fipscode,
                       entrant,
                       [entrant],
                       DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                       WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                       WeightVec([0.1, 0.1, 0.1]), #TODO: needs to be fixed with logitest
                       neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0),
                       Array{Int64, 1}())
      push!(el.config, entr) # need to create a new record for this hospital in the market
      # need to add it to the dictionary too:
      el.collection[newfid] = entr
      for elm in el.config
        NeighborAppend(elm, entr)
        NeighborAppend(entr, elm)
      end

    end
    for elm in el.config
    # This does the actual sampling process - Takes the hospital in elm, selects the corresponding choices, then samples according to the probs.
  #    println(sample( ChoicesAvailable(elm), elm.chprobability ))
      newchoice = LevelFunction(elm, sample( ChoicesAvailable(elm), elm.chprobability ))
      elm.level = newchoice
      push!(elm.levelhistory, newchoice)
    #  println(elm.levelhistory)
    end
  end
end
for el in Texas.ms
  MarketCleaner(el) # Remove Entrants from Market Record.
end








###
