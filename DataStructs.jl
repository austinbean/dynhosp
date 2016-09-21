
# Store the hospitals in data structures to reduce the complexity of writing out all of the values
# to matrices and keeping track of all of the indices, etc.

include("/Users/austinbean/Desktop/dynhosp/Reboot.jl")

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
  probhistory::Array{Float64,1}
  # The logitest function takes the following:
  # logitest((0,0), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
  neigh::neighbors
  hood::Array{Int64, 1}
  bedcount::Float64
  perturbed::Bool
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
# Texas = EntireState(Array{hospital,1}(), Dict{Int64, Market}(), Dict{Int64, hospital}())
# See below for dictionary comprehension.

fips = unique(data[:,78])
# This adds a list of the markets covered to the whole state iterable.
function MakeIt(Tex::EntireState, fip::Vector)
  for el in fip
    if el != 0
      el = eval(parse("m$el = Market( Array{hospital,1}(), Dict{Int64, hospital}() ,$el)"))
      push!(Tex.ms, el)
    end
  end
  Tex.mkts = [ m.fipscode => m for m in Tex.ms]
end

#MakeIt(Texas, fips)
# Write the market dictionary out as a comprehension:
# Texas.mkts = [ m.fipscode => m for m in Texas.ms]

# fid - col 74, lat - col 94, long - col 95, name - col 82, fipscode - col 78, act_int - col 79, act_solo - col 80
data05 = data[(data[:,75].==2005), :] ;

function TXSetup(Tex::EntireState, data::Matrix; lev105loc = 97, lev205loc = 98, lev305loc = 99, lev1515loc = 101, lev2515loc = 102, lev3515loc = 103, lev11525loc = 105, lev21525loc = 106, lev31525loc = 107)
  for i = 1:size(data,1)
    fips = data[i, 78]
    if fips != 0
      level = 0
      if (data[i, 79] == 1)&(data[i,80]==0)
        level = 3
      elseif (data[i, 79] == 0)&(data[i,80]==1)
        level = 2
      else
        level = 1
      end
      push!(Tex.mkts[fips].config,
      hospital( data[i, 74],
                data[i,94],
                data[i, 95],
                data[i, 82],
                fips,
                level,
                [level],
                DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                WeightVec([data[i,19], data[i,37], data[i,55], data[i, 73]]),
                Array{Float64,1}(),
                neighbors(data[i, lev105loc], data[i,lev205loc ], data[i,lev305loc ], data[i,lev1515loc ], data[i,lev2515loc ], data[i, lev3515loc], data[i,lev11525loc ], data[i,lev21525loc ], data[i,lev31525loc]  ),
                Array{Int64,1}(),
                  0    , # need the beds here - currently NOT in this data.
                false ) )
    end
    # push all hospital fid/ fips pairs into the directory.
    Tex.fipsdirectory[data[i, 74]] = fips # now for the whole state I can immediately figure out which market a hospital is in.
  end
  return Tex
end

#Texas = TXSetup(Texas, data05);

# Expand the market dictionaries so that they are filled with the hospitals
function ExpandDict(Tex::EntireState)
  for el in Tex.ms
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
end

#ExpandDict(Texas)

function MakeNew(fi::Vector, dat::Matrix)
  Texas = EntireState(Array{hospital,1}(), Dict{Int64, Market}(), Dict{Int64, hospital}())
  MakeIt(Texas, fi)
  TXSetup(Texas, dat)
  ExpandDict(Texas)
  return Texas
end

Texas = MakeNew(fips, data05);


function MarketPrint(mkt::Market)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  println(mkt.fipscode)
  for el in mkt.config
    println("*******************")
    println(el.name)
    println(el.neigh)
    println(el.hood)
  end
end

function NeighborsPrint(mkt::Market)
  for el in mkt.config
    println(el.name, " ", el.neigh)
  end
end

function FacPrint(hosp::hospital)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  println(hosp.name)
  println(hosp.fid)
  println("Fips: ", hosp.fipscode)
  println("Level: ", hosp.level)
  println("Neighbors: ", hosp.hood)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
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




function MktSize(n::neighbors)
  sum1 = n.level105 + n.level1515 + n.level11525
  sum2 = n.level205 + n.level2515 + n.level21525
  sum3 = n.level305 + n.level3515 + n.level31525
  return sum1, sum2, sum3
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
    if (dist < 25)&(entrant.level != -999)
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
  takes two hospital records, computes the distance between them and subtracts 1 from the relevant record in the neighborhood type.
  It removes the record of entrant FROM the record of elm.
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
    deleteat!(elm.hood, findin(elm.hood, entrant.fid))
  end
end

function NeighborClean(state::EntireState)
  # This will set every value in the neighbors category of every hospital in the state to zero.
  for mkt in state.ms
    for hosp in mkt.config
      hosp.neigh = neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0)
      hosp.hood = Array{Int64,1}()
    end
  end
end


function NeighborFix(state::EntireState)
  # For every hospital in the state, append all other hospitals within 25 miles, ignoring county boundaries
  for mkt1 in state.ms
    for mkt2 in state.ms
      if mkt1.fipscode != mkt2.fipscode
        for hos1 in mkt1.config
          for hos2 in mkt2.config
            if distance(hos1.lat, hos1.long, hos2.lat, hos2.long) < 25
              NeighborAppend(hos1, hos2)
              NeighborAppend(hos2, hos1)
            end
          end
        end
      end
    end
  end
end

#NB: The function below will fix the neighbors while respecting the county boundaries, unlike the above.

function StrictCountyNeighborFix(state::EntireState)
  # For every hospital in the state, append all other hospitals within 25 miles
  for mkt1 in state.ms
    for hos1 in mkt1.config
      for hos2 in mkt1.config
        if hos1.fid != hos2.fid
          if distance(hos1.lat, hos1.long, hos2.lat, hos2.long) < 25
            NeighborAppend(hos1, hos2)
            NeighborAppend(hos2, hos1)
          end
        end
      end
    end
  end
end

function HospFindFirst(mkt::Market, hosp::hospital)
  # looks for a hospital given by hosp in the market mkt, by searching for the fid and returning the index.
  found = 0
  for el in 1:size(mkt.config,1)
    if mkt.config[el].fid == hosp.fid
      found = el
    end
  end
  return found
end

function FidFindFirst(mkt::Market, fid::Int64)
  # looks for a fid in the market, then returns the index of the fid.
  found = 0
  for el in 1:size(mkt.config,1)
    if mkt.config[el].fid == fid
      found = el
    end
  end
  return found
end

function MarketCleaner(mkt::Market)
  # Should take as an argument a whole market and then remove the records of any entrants.
  entlist = Array{Int64,1}()
  for el in mkt.config
    if el.fid < 0
      push!(entlist, el.fid)
    end
  end
  for hosps in mkt.config
    for exfid in entlist
      if in(exfid, hosps.hood)
        NeighborRemove(hosps, mkt.config[FidFindFirst(mkt, exfid)] ) #TODO: This is not getting the distance change right.
      end
    end
  end
  for el in entlist
  #  println(el)
    # Remove the hospital from the market array
    deleteat!(mkt.config, FidFindFirst(mkt, el)) # NB: HospFindFirst takes *market* as argument, but deleteat! takes *array*, i.e, market.config
    # Remove the hospital from the market dictionary
    pop!(mkt.collection, el)
  end
end


function HospUpdate(hosp::hospital, choice::Int)
  levl = (-1, -1)
 if (hosp.level!=choice)
   if choice != -999
     if choice == 1
       levl = (0,0)
     elseif choice == 2
       levl = (1,0)
     elseif choice == 3
       levl = (0,1)
     end
     levels = MktSize(hosp.neigh)
     prs = logitest(levl, levels[1], levels[2], levels[3], [hosp.neigh.level105; hosp.neigh.level205; hosp.neigh.level305; hosp.neigh.level1515; hosp.neigh.level2515; hosp.neigh.level3515; hosp.neigh.level11525; hosp.neigh.level21525; hosp.neigh.level31525 ] )
  #   println(prs)
     return WeightVec(vec(prs))
   else # choice = -999
     return WeightVec([1.0]) #TODO: one option, no choices ??  Might need four options [1.0 1.0 1.0 1.0]
   end
 else
   return hosp.chprobability
  end
end


function NewSim(T::Int, Tex::EntireState; entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002] )
  for i = 1:T
    for el in Tex.ms
      entrant = sample(entrants, WeightVec(entryprobs))
      if entrant != 0
        println("Entry! ", entrant , " FIPS ", el.fipscode)
        entloc = NewEntrantLocation(el) # called on the market
        newfid = -floor(rand()*1e6)
        entr = hospital( newfid, entloc[1], entloc[2], " Entrant $newfid ", el.fipscode, entrant, [entrant],
                         DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                         WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                         WeightVec([0.1, 0.1, 0.1, 0.1]), Array{Float64,1}(), neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0), Array{Int64, 1}(), 0, false)
        push!(el.config, entr) # need to create a new record for this hospital in the market
        # need to add it to the dictionary too:
        el.collection[newfid] = entr
        for elm in el.config
          NeighborAppend(elm, entr)
          NeighborAppend(entr, elm)
        end
        # TODO: Call Hospital Update here once the neighbors exist
         HospUpdate(entr, entrant) #entrant is the level
      end
      for elm in el.config
      #  println("Facility ", elm.fid, "   ", elm.fipscode)
        #TODO: There may be a very infrequent error that pops up connected to HospUpdate.
        action = sample( ChoicesAvailable(elm), elm.chprobability )                            # Take the action
        push!( elm.probhistory ,elm.chprobability[ findin(ChoicesAvailable(elm), action)[1] ]) # Record the prob with which the action was taken.
        newchoice = LevelFunction(elm, action)                                                 # What is the new level?
        elm.chprobability = HospUpdate(elm, newchoice)                                         # What are the new probabilities, given the new level?
        elm.level = newchoice                                                                  # Set the level to be the new choice.
        push!(elm.levelhistory, newchoice)
      end
    end
    NeighborClean(Tex)
    NeighborFix(Tex)
  end
end

NewSim(50, Texas)

for el in Texas.ms
  MarketCleaner(el) # Remove Entrants from Market Record.
end



type patientcount
 count385::Int64
 count386::Int64
 count387::Int64
 count388::Int64
 count389::Int64
 count390::Int64
 count391::Int64
end


type coefficients
  distance::Float64
  distsq::Float64
  inten::Float64
  inter::Float64
  distbed::Float64
  closest::Float64
  # can add extras
end

type zip
 code::Int64
 phr::Int64 # may have coefficients differing by PHR
 facilities::Dict{Int64, hospital}
 fes::Dict{Int64, Float64} # keep a dict of hospital FE's at the zip around
 pdetutils::Dict{Int64, Float64} # keep the deterministic utilities
 mdetutils::Dict{Int64, Float64} # the same for medicare patients.
 lat::Float64
 long::Float64
 pcoeffs::coefficients
 mcoeffs::coefficients
 ppatients::patientcount
 mpatients::patientcount
end

type patientcollection
 zips::Dict{Int64, zip}
end

function CreateZips(zipcodes::Array, ch::Array, Tex::EntireState; phrloc = 103)
  ppatients = patientcollection( Dict{Int64, zip}() )
  unfound = Array{Int64,1}()
  for el in zipcodes
    ppatients.zips[el] = zip(el, 0, Dict{Int64,hospital}(), Dict{Int64,Float64}(), # zipcode, public health region, facilities, hospital FE's.
                             Dict{Int64,Float64}(), Dict{Int64, Float64}(), # private det utilities, medicaid det utilities.
                             0.0, 0.0, coefficients(privatedistance_c, privatedistsq_c, privateneoint_c, privatesoloint_c, privatedistbed_c, privateclosest_c),
                             coefficients(medicaiddistance_c, medicaiddistsq_c, medicaidneoint_c, medicaidsoloint_c, medicaiddistbed_c, medicaidclosest_c),  #lat, long, private coefficients, medicaid coefficients
                             patientcount(0,0,0,0,0,0,0), patientcount(0,0,0,0,0,0,0)) # private and medicaid patients
  end
  for i = 1:size(ch, 1) #rows
    ppatients.zips[ch[i,1]].lat = ch[i,7]
    ppatients.zips[ch[i,1]].long = ch[i,8]
    for j = 11:17:size(ch,2)
      try
        Tex.fipsdirectory[ch[i,j]]
        # NB: choices[i,11] → FID, Tex.fipsdirectory: FID → FIPSCODE.
        # NB: Tex.fipsdirectory[ ch[i, 11]]: FID → Market,
        # NB: Tex.fipsdirectory[ ch[i, 11]].collection[ ch[i, 11]]: FID → Hospital Record.
        fipscode = Tex.fipsdirectory[ch[i,j]]
        ppatients.zips[ch[i,1]].facilities[ch[i,j]] = Tex.mkts[fipscode].collection[ch[i,j]]
        Tex.mkts[ fipscode].collection[ ch[i, j]].bedcount = ch[i,j+3]
      catch y
        if isa(y, KeyError)
          push!(unfound,ch[i, j])
        end
      end
    end
  end
  return ppatients, unfound
end

patients, unf = CreateZips(zips, choices, Texas);


function PrintZip(zi::zip)
  # Prints the fid and the name of the facilities attached to the zips.
  for el in keys(zi.facilities)
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
    println(el, "  ", zi.facilities[el].name)
  end
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
    println("385 ", zi.ppatients.count385, " ", zi.mpatients.count385)
    println("386 ", zi.ppatients.count386, " ", zi.mpatients.count386)
    println("387 ", zi.ppatients.count387, " ", zi.mpatients.count387)
    println("388 ", zi.ppatients.count388, " ", zi.mpatients.count388)
    println("389 ", zi.ppatients.count389, " ", zi.mpatients.count389)
    println("390 ", zi.ppatients.count390, " ", zi.mpatients.count390)
    println("391 ", zi.ppatients.count391, " ", zi.mpatients.count391)
end

function FillPPatients(pats::patientcollection, imported::Matrix; ziploc = 101, drgloc = 104)
  notfound = Array{Int64,1}()
  for row in 1:size(imported, 1)
    if imported[row, drgloc ] == 385
      pats.zips[imported[row, ziploc]].ppatients.count385 += 1;
    elseif imported[row, drgloc ] == 386
      pats.zips[imported[row, ziploc]].ppatients.count386 += 1;
    elseif imported[row, drgloc ] == 387
      pats.zips[imported[row, ziploc]].ppatients.count387 += 1;
    elseif imported[row, drgloc ] == 388
      pats.zips[imported[row, ziploc]].ppatients.count388 += 1;
    elseif imported[row, drgloc ] == 389
      pats.zips[imported[row, ziploc]].ppatients.count389 += 1;
    elseif imported[row, drgloc ] == 390
      pats.zips[imported[row, ziploc]].ppatients.count390 += 1;
    elseif imported[row, drgloc ] == 391
      pats.zips[imported[row, ziploc]].ppatients.count391 += 1;
    else # not found?
        push!(notfound, pats.zips[imported[row, ziploc]].code);
  #     println( imported[row, drgloc])
    end
  end
  return pats;
end

function FillMPatients(pats::patientcollection, imported::Matrix; ziploc = 101, drgloc = 104)
  notfound = Array{Int64,1}()
  for row in 1:size(imported, 1)
    if imported[row, drgloc ] == 385
      pats.zips[imported[row, ziploc]].mpatients.count385 += 1;
    elseif imported[row, drgloc ] == 386
      pats.zips[imported[row, ziploc]].mpatients.count386 += 1;
    elseif imported[row, drgloc ] == 387
      pats.zips[imported[row, ziploc]].mpatients.count387 += 1;
    elseif imported[row, drgloc ] == 388
      pats.zips[imported[row, ziploc]].mpatients.count388 += 1;
    elseif imported[row, drgloc ] == 389
      pats.zips[imported[row, ziploc]].mpatients.count389 += 1;
    elseif imported[row, drgloc ] == 390
      pats.zips[imported[row, ziploc]].mpatients.count390 += 1;
    elseif imported[row, drgloc ] == 391
      pats.zips[imported[row, ziploc]].mpatients.count391 += 1;
    else # not found?
        push!(notfound, pats.zips[imported[row, ziploc]].code);
  #     println( imported[row, drgloc])
    end
  end
  return pats;
end

function FillPatients(pats::patientcollection, private::Matrix, medicaid::Matrix)
  pats = FillPPatients(pats, private)
  pats = FillMPatients(pats, medicaid)
end

patients = FillPatients(patients, pinsured, pmedicaid);

function ComputeDetUtil(zipc::zip, fid::Int64, p_or_m::Bool)
  # Computes the deterministic component of utility for each hospital.
  dist = distance(zipc.facilities[fid].lat, zipc.facilities[fid].long, zipc.lat, zipc.long)
  if p_or_m #if TRUE private
    if zipc.facilities[fid].level == 1
      return zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 2
      return zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.inter+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 3
      return zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.inten+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    else  # =-999
      return - 99 # can't choose a facility which has exited - set det utility very low.
    end
  else
    if zipc.facilities[fid].level == 1
      return zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 2
      return zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.inter+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 3
      return zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.inten+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    else  # =-999
      return - 99
    end
  end
end

function UpdateDeterministic(collect::patientcollection)
  for el in keys(collect.zips) #iterates over zips
    for fid in keys(collect.zips[el].facilities) # iterates over dict of facilities within zip.
      if (size(collect.zips[el].facilities[fid].levelhistory,1)==1)||(collect.zips[el].facilities[fid].level!=collect.zips[el].facilities[fid].levelhistory[end]) #state has changed OR this is the first period.
        collect.zips[el].mdetutils[fid] = ComputeDetUtil(collect.zips[el], fid, false)
        collect.zips[el].pdetutils[fid] = ComputeDetUtil(collect.zips[el], fid, true)
      end
    end
  end
end

function GenChoices(collect::zip; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  for el in keys(zip.facilities) # dict has field count - has element number. 
      utils = hcat([ [k1,collect.zips[el].pdetutils[k1]] for k1 in keys(collect.zips[el].pdetutils)]...)
      temparr = zeros(size(utils, 2))
      for k = 1:collect.zips[el].ppatients.count385
        utils[1,indmax(utils + rand!(d, temparr))] += 1
      end
      for k = 1:collect.zips[el].ppatients.count386
        utils[1,indmax(utils + rand!(d, temparr))] += 1
      end
      for k = 1:collect.zips[el].ppatients.count387
        utils[1,indmax(utils + rand!(d, temparr))] += 1
      end
      for i=1:collect.zips[el].ppatients.count388
        utils[1,indmax(utils + rand!(d, temparr))] += 1
      end
      for i = 1:collect.zips[el].ppatients.count389
        utils[1,indmax(utils + rand!(d, temparr))] += 1
      end
      for i=1:collect.zips[el].ppatients.count390
        utils[1,indmax(utils + rand!(d, temparr))] += 1
      end
      for i = 1:collect.zips[el].ppatients.count391
        utils[1,indmax(utils + rand!(d, temparr))] += 1
      end
  end
end

###
