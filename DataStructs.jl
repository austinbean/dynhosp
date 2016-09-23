
# Store the hospitals in data structures to reduce the complexity of writing out all of the values
# to matrices and keeping track of all of the indices, etc.

include("/Users/austinbean/Desktop/dynhosp/Reboot.jl")

          ### NB: Supply-side Data Structures. ###
type WTP
  w385::Array{Float64, 1}
  w386::Array{Float64, 1}
  w387::Array{Float64, 1}
  w388::Array{Float64, 1}
  w389::Array{Float64, 1}
  w390::Array{Float64, 1}
  w391::Array{Float64, 1}
end

type DemandHistory
  demand385::Array{Int64, 1}
  demand386::Array{Int64, 1}
  demand387::Array{Int64, 1}
  demand388::Array{Int64, 1}
  demand389::Array{Int64, 1}
  demand390::Array{Int64, 1}
  demand391::Array{Int64, 1}
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
  fipscode::Int64
  level::Int64
  levelhistory::Vector{Int64}
  pdemandhist::DemandHistory # separate histories for Private and Medicaid patients.
  mdemandhist::DemandHistory
  wtphist::WTP
  chprobability::WeightVec
  probhistory::Array{Float64,1}
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

        #### NB: Demand-side Data Structures ######

type patientcount
 count385::Int64
 count386::Int64
 count387::Int64
 count388::Int64
 count389::Int64
 count390::Int64
 count391::Int64
end

import Base.+
function +(x::patientcount, y::patientcount)
  return patientcount(x.count385 + y.count385, x.count386 + y.count386, x.count387 + y.count387, x.count388 + y.count388, x.count389 + y.count389, x.count390 + y.count390, x.count391 + y.count391)
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



        ##### NB: Supply-side Data Creation Functions ######

function MakeIt(Tex::EntireState, fip::Vector)
  # Perhaps poor practice to use Eval in this way, but generates markets named m*fipscode* for any fipscode in the vector fip.
  for el in fip
    if el != 0
      el = eval(parse("m$el = Market( Array{hospital,1}(), Dict{Int64, hospital}() ,$el)"))
      push!(Tex.ms, el)
    end
  end
  Tex.mkts = [ m.fipscode => m for m in Tex.ms]
end


function TXSetup(Tex::EntireState, data::Matrix; lev105loc = 97, lev205loc = 98, lev305loc = 99, lev1515loc = 101, lev2515loc = 102, lev3515loc = 103, lev11525loc = 105, lev21525loc = 106, lev31525loc = 107)
  # Takes an entire state and adds data from the imported choices returns a record with
  # fipscodes containing hospitals with mostly empty field values.
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
                DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                WeightVec([data[i,19], data[i,37], data[i,55], data[i, 73]]),
                Array{Float64,1}(),
                neighbors(data[i, lev105loc], data[i,lev205loc ], data[i,lev305loc ], data[i,lev1515loc ], data[i,lev2515loc ], data[i, lev3515loc], data[i,lev11525loc ], data[i,lev21525loc ], data[i,lev31525loc]  ),
                Array{Int64,1}(),
                  0    , # beds added later.
                false ) )
    end
    # push all hospital fid/ fips pairs into the directory.
    Tex.fipsdirectory[data[i, 74]] = fips # now for the whole state I can immediately figure out which market a hospital is in.
  end
  return Tex
end


function ExpandDict(Tex::EntireState)
  # Expand the market dictionaries so that they are filled with the hospitals
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


function MakeNew(fi::Vector, dat::Matrix)
  # Call this and the whole state with all markets should be created.
  Texas = EntireState(Array{hospital,1}(), Dict{Int64, Market}(), Dict{Int64, hospital}())
  MakeIt(Texas, fi)
  TXSetup(Texas, dat)
  ExpandDict(Texas)
  return Texas
end

# Data - should be moved to Reboot.jl eventually.
fips = unique(data[:,78])
data05 = data[(data[:,75].==2005), :] ;

#NB: Creates hospital datastructure
Texas = MakeNew(fips, data05);

      #### NB:  Supply-side  Printing Utilities to Display Simulation Outcomes
function MarketPrint(mkt::Market)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  println(mkt.fipscode)
  for el in mkt.config
    println("*******************")
    println(el.name)
    println(el.neigh)
    println(el.hood)
    println(el.chprobability)
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

  ### NB: Substantive Supply-side Functions.

function NewEntrantLocation(mkt::Market)
  # Takes the market, takes the mean location of all hospitals, adds normal noise to it.  ≈ 6 miles perturbation from mean.
  meanlat = 0
  meanlong = 0
  for el in mkt.config # over hospitals
    meanlat += el.lat
    meanlong += el.long
  end
  return [meanlat/size(mkt.config, 1) + rand(Normal(0, 0.1), 1)[1], meanlong/size(mkt.config, 1) + rand(Normal(0, 0.1), 1)[1]]
end


function MktSize(n::neighbors)
  # takes a set of neighbors and returns the sum of levels 1, 2, 3 at the various distances.
  sum1 = n.level105 + n.level1515 + n.level11525
  sum2 = n.level205 + n.level2515 + n.level21525
  sum3 = n.level305 + n.level3515 + n.level31525
  return sum1, sum2, sum3
end

function ChoicesAvailable(h::hospital)
  # Takes a hospital, returns the choices available at that level as vectors.
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
  Takes two hospital records, computes the distance between them and adds a 1 to the relevant record in the neighborhood type.
  Appends it to the hood of elm, which is a list of fids.  So this adds to both elm.neigh and elm.hood.
  It is not symmetric - it appends entrant to elm, not vice versa.
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
  It removes the record of entrant FROM the record of elm.  Also not symmetric - removes entrant from elm's records, not the reverse.
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


function StrictCountyNeighborFix(state::EntireState)
  # For every hospital in the state, append all other hospitals within 25 miles AND in the same county.
  # More restrictive than NeighborFix
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
  # Takes a whole market and then removes the records of any entrants.
  entlist = Array{Int64,1}()
  for el in mkt.config
    if el.fid < 0 # all entrants are tagged with negative fids.
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
    # Remove the hospital from the market array
    deleteat!(mkt.config, FidFindFirst(mkt, el)) # NB: HospFindFirst takes *market* as argument, but deleteat! takes *array*, i.e, market.config
    # Remove the hospital from the market dictionary
    pop!(mkt.collection, el)
  end
end


function HospUpdate(hosp::hospital, choice::Int; update = false)
  # Takes a hospital record and updates the probabilities of the choices.
  levl = (-1, -1)
 if (hosp.level!=choice)|update # want to be able to force this to rerun when the data is cleaned again.
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

      ##### NB: Demand-side Data Structure Creation.


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

function FillPPatients(pats::patientcollection, imported::Matrix; ziploc = 101, drgloc = 104)
  # Takes the imported matrix of *privately-insured* patients and records the number at each DRG 385-391 in each zip record.
  # There is a separate function for the Medicaid patients.
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
    end
  end
  return pats;
end

function FillMPatients(pats::patientcollection, imported::Matrix; ziploc = 101, drgloc = 104)
  # Takes the imported matrix of *Medicaid* patients and records the number at each DRG 385-391 in each zip record.
  # There is a separate function for the privately-insured patients.
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
    end
  end
  return pats;
end

function FillPatients(pats::patientcollection, private::Matrix, medicaid::Matrix)
  # Adds the privately insured and medicaid patients to the zip records.
  pats = FillPPatients(pats, private)
  pats = FillMPatients(pats, medicaid)
end

patients = FillPatients(patients, pinsured, pmedicaid);


    ### NB: Zip code record printing utility.

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

  ### NB: Substantive Demand-side Functions.

  #TODO: Modify to add the outside option!

function ComputeDetUtil(zipc::zip, fid::Int64, p_or_m::Bool)
  # Computes the deterministic component of utility for each hospital in the zip "zipc".
  # Maps exited facilites to have deterministic utility -999
  # Works on private and medicaid patients by setting p_or_m to true or false, respectively.
  # Has been written to accomodate hospital FE's when available.
  dist = distance(zipc.facilities[fid].lat, zipc.facilities[fid].long, zipc.lat, zipc.long)
  if p_or_m #if TRUE private
    if zipc.facilities[fid].level == 1
      return zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 2
      return zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.inter+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 3
      return zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.inten+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    else  # =-999
      return -999 # can't choose a facility which has exited - set det utility very low.
    end
  else
    if zipc.facilities[fid].level == 1
      return zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 2
      return zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.inter+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    elseif zipc.facilities[fid].level == 3
      return zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.inten+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)+zipc.pcoeffs.closest*(0) #+ zipc.fes[fid]
    else  # =-999
      return -999
    end
  end
end

function WhichZips(pats::patientcollection, fid::Int64)
  # Takes a patientcollection and tells me which zips have the hospital fid
  for zi in keys(pats.zips)
    try
      pats.zips[zi].facilities[fid]
      println(zi)
    catch y
      if isa(y, KeyError)
        #not found.  Whatever.
      end
    end
  end
end



function CalcWTP(zipc::zip)
  # Takes the deterministic component of utility for the privately insured patients and returns a WTP measure.
  outp = [j => 0.0 for j in keys(zipc.pdetutils)]
  interim = 0.0
  for el in keys(zipc.pdetutils)
    #NB: computed utility of exited firm will be zero by ComputeDetUtil assigning it to -999, so exp(-999) = 0
    interim +=  (outp[el] = exp(zipc.pdetutils[el]) ) #NB: This is a nice trick - simultaneously assigning and adding.
  end
  return [ j => outp[j]/interim for j in keys(outp)]
end


#TODO: here we need to make sure that there is always an outside option which is zeros.
#TODO: Also - these need to vary by DRG.  That requires: detutil by the DRG at the zipcode level.

function WTPMap(pats::patientcollection, Tex::EntireState)
  # Takes a patient collection and an entire state and returns a dict{fid, WTP}
  # computed by calling CalcWTP.  Right now it ignores Inf and NaN.
  outp = [ j=> 0.0 for j in keys(Tex.fipsdirectory) ]
  for zipc in keys(pats.zips)
    vals = CalcWTP(pats.zips[zipc])
    for el in keys(vals) # What to do about key errors?  there will be some.
      try
        outp[el]
        if (vals[el]!=1)&!isnan(vals[el])
          outp[el]+= (1/(1-vals[el]))
        end
      catch y # here the issue was that this was catching an "InexactError" but there was no test for it.
        if isa(y, KeyError)
          #println(el) #the facility is missing.  Could be an entrant.
        end
      end
    end
  end
  return outp # gives a dict{fid, WTP} back
end

function WriteWTP(reslt::Dict{Int64, Float64}, Tex::EntireState)
  # Takes a dict of {fid, WTP} and writes it out by DRG.
  for els in keys(reslt)
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w385, reslt[els])
  #  push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w386, reslt[els])
  #  push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w387, reslt[els])
  #  push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w388, reslt[els])
  #  push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w389, reslt[els])
  #  push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w390, reslt[els])
  #  push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w391, reslt[els])
  end
end

function UpdateDeterministic(collect::patientcollection)
  # Computes the deterministic component of the utility - updates every firm every time it is called.
  for el in keys(collect.zips) #iterates over zips
    for fid in keys(collect.zips[el].facilities) # iterates over dict of facilities within zip.
      collect.zips[el].mdetutils[fid] = ComputeDetUtil(collect.zips[el], fid, false)
      collect.zips[el].pdetutils[fid] = ComputeDetUtil(collect.zips[el], fid, true)
    end
  end
end

UpdateDeterministic(patients)

function GenPChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  # The patient choice is max \bar{U} + ϵ, but we have \bar{U} from Compute Det Util and we know how many patients there are in
  # the privately insured category from FillPPatients.  This returns a dict of fids and patient counts, where patient counts are
  # generated by repeatedly finding max i = 1, ..., N \bar{U}_i + ϵ_i.  Note the corresponding GenMChoices below.
  outp = [ j => patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory)] # output is a {FID, patientcount} dictionary.
  for zipcode in keys(pats.zips)
    if pats.zips[zipcode].pdetutils.count > 0
      utils = hcat([ [k1,pats.zips[zipcode].pdetutils[k1]] for k1 in keys(pats.zips[zipcode].pdetutils)]...)
      temparr = zeros(size(utils, 2))
      for k = 1:pats.zips[zipcode].ppatients.count385
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count385 += 1
      end
      for k = 1:pats.zips[zipcode].ppatients.count386
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count386 += 1
      end
      for k = 1:pats.zips[zipcode].ppatients.count387
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count387 += 1
      end
      for i=1:pats.zips[zipcode].ppatients.count388
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count388 += 1
      end
      for i = 1:pats.zips[zipcode].ppatients.count389
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count389 += 1
      end
      for i=1:pats.zips[zipcode].ppatients.count390
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count390 += 1
      end
      for i = 1:pats.zips[zipcode].ppatients.count391
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count391 += 1
      end
    end
  end
  return outp
end


function GenMChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  # The patient choice is max \bar{U} + ϵ, but we have \bar{U} from Compute Det Util and we know how many patients there are in
  # the Medicaid category from FillMPatients.  This returns a dict of fids and patient counts, where patient counts are
  # generated by repeatedly finding max i = 1, ..., N \bar{U}_i + ϵ_i.  Note the corresponding GenPChoices above.
  outp = [ j => patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory)] # output is a {FID, patientcount} dictionary.
  for zipcode in keys(pats.zips)
    if pats.zips[zipcode].mdetutils.count > 0
      utils = hcat([ [k1,pats.zips[zipcode].mdetutils[k1]] for k1 in keys(pats.zips[zipcode].mdetutils)]...)
      temparr = zeros(size(utils, 2))
      for k = 1:pats.zips[zipcode].mpatients.count385
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count385 += 1
      end
      for k = 1:pats.zips[zipcode].mpatients.count386
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count386 += 1
      end
      for k = 1:pats.zips[zipcode].mpatients.count387
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count387 += 1
      end
      for i=1:pats.zips[zipcode].mpatients.count388
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count388 += 1
      end
      for i = 1:pats.zips[zipcode].mpatients.count389
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count389 += 1
      end
      for i=1:pats.zips[zipcode].mpatients.count390
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count390 += 1
      end
      for i = 1:pats.zips[zipcode].mpatients.count391
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr)')]].count391 += 1
      end
    end
  end
  return outp
end

function PHistoryAdd(hos::hospital, cnt::patientcount)
  # Maps patientcount to the private demand history
  push!(hos.pdemandhist.demand385, cnt.count385)
  push!(hos.pdemandhist.demand386, cnt.count386)
  push!(hos.pdemandhist.demand387, cnt.count387)
  push!(hos.pdemandhist.demand388, cnt.count388)
  push!(hos.pdemandhist.demand389, cnt.count389)
  push!(hos.pdemandhist.demand390, cnt.count390)
  push!(hos.pdemandhist.demand391, cnt.count391)
end

function MHistoryAdd(hos::hospital, cnt::patientcount)
  # Maps patientcount to the Medicaid demand history.
  push!(hos.mdemandhist.demand385, cnt.count385)
  push!(hos.mdemandhist.demand386, cnt.count386)
  push!(hos.mdemandhist.demand387, cnt.count387)
  push!(hos.mdemandhist.demand388, cnt.count388)
  push!(hos.mdemandhist.demand389, cnt.count389)
  push!(hos.mdemandhist.demand390, cnt.count390)
  push!(hos.mdemandhist.demand391, cnt.count391)
end

function PDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)
  for el in keys(patd)
    PHistoryAdd(Tex.mkts[Tex.fipsdirectory[el]].collection[el], patd[el])
  end
end

function MDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)
  for el in keys(patd)
    MHistoryAdd(Tex.mkts[Tex.fipsdirectory[el]].collection[el], patd[el])
  end
end



function HospitalClean(hos::hospital)
  # resets the hospital to the initial state after a run of the simulation.
  # some fields don't change: fid, lat, long, name, fips, bedcount,
  # eventually perturbed will probably change.
  hos.level = hos.levelhistory[1]                  #Reset level history to initial value
  hos.levelhistory = [hos.levelhistory[1]]           #Set record length back to zero.
  hos.mdemandhist = DemandHistory( Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}()) #empty demand history
  hos.pdemandhist = DemandHistory( Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}()) #empty demand history
  hos.wtphist = WTP(Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ) #empty WTP
  hos.chprobability = WeightVec([0])               #TODO: this one is complicated!
  hos.probhistory = Array{Float64,1}()             #Empty history of choices.
  hos.neigh = neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0) #TODO: reset this in Restore()
  hos.hood = Array{Int64,1}()                      #TODO: reset in Restore()
  hos.perturbed = false                            #For now this is always false.
end

function Restore(Tex::EntireState)
  # This function needs to set all hospital states back to zero.
  # Then it re-computes the set of neighbors using NeighborFix, fixing hosp.neigh and hosp.hood.
  # Finally it recomputes the initial choice probabilities.
  #TODO - there are entrants in this group too.
  for mkt in Tex.ms
    for hos in mkt.config
      HospitalClean(hos)
    end
  end
#TODO: this isn't quite working yet.  But it isn't crucial at the moment.

  NeighborFix(Tex) # Restores all neighbors to both hosp.neigh and hosp.hood.
  for mkt in Tex.ms
    for hos in mkt.config
      HospUpdate(hos, hos.level; update = true) # HospUpdate should now fix these.
    end
  end
  return Tex
end



    ### NB: The business of the simulation.



function NewSim(T::Int, Tex::EntireState, pats::patientcollection; entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002] )
  for i = 1:T
    println(i)
    WriteWTP(WTPMap(pats, Tex), Tex)
    PDemandMap(GenPChoices(pats, Tex), Tex)
    MDemandMap(GenMChoices(pats, Tex), Tex)
    for el in Tex.ms
      entrant = sample(entrants, WeightVec(entryprobs))
      if entrant != 0
      #  println("Entry! ", entrant , " FIPS ", el.fipscode)
        entloc = NewEntrantLocation(el) # called on the market
        newfid = -floor(rand()*1e6) # all entrant fids negative to facilitate their removal later.
        entr = hospital( newfid, entloc[1], entloc[2], " Entrant $newfid ", el.fipscode, entrant, [entrant],
                         DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
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
         HospUpdate(entr, entrant) #entrant is the level
      end
      for elm in el.config
        action = sample( ChoicesAvailable(elm), elm.chprobability )                            # Take the action
        push!( elm.probhistory ,elm.chprobability[ findin(ChoicesAvailable(elm), action)[1] ]) # Record the prob with which the action was taken.
        newchoice = LevelFunction(elm, action)                                                 # What is the new level?
        elm.chprobability = HospUpdate(elm, newchoice)                                         # What are the new probabilities, given the new level?
        elm.level = newchoice                                                                  # Set the level to be the new choice.
        push!(elm.levelhistory, newchoice)
      end
    end
    # This updates after all of the facilities have been changed.
      UpdateDeterministic(pats)                                                                # Updates deterministic component of utility
      for zipc in keys(pats.zips)
        pdetutils = CalcWTP(pats.zips[zipc])                                                   # Calculates WTP from deterministic utility, now updated.
      end
    #TODO: Figure out the cleanup things, but they are turned off for now.  Alternatively - remake every time in the larger sim.
  #  NeighborClean(Tex)
  #  NeighborFix(Tex)
  end
  return Tex
end

NewSim(10, Texas, patients)

for el in Texas.ms
  MarketCleaner(el) # Remove Entrants from Market Record.
end


function OuterSim(MCcount::Int)
  MakeIt() #recreate the state
  #TODO: call the function to create the individual records.
  for j = 1:MCcount
    NTex = NewSim(50, Texas, patients) # generates the sim results.
    ResultsOut(NTex)
  end

end

function ResultsOut(Tex::EntireState)
  for el in keys(Tex.fipsdirectory) # Now this is all of the hospitals.

  end
end





### TODO: What's next?

#=

- Outer function for MC sims.
- Map Results to payoff function.  This first.
- Fix WTP.


=#




###
