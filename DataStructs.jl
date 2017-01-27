#=

Works on this version if you change the dictionary generation syntax:
Julia Version 0.4.5
Commit 2ac304d (2016-03-18 00:58 UTC)
Platform Info:
  System: Darwin (x86_64-apple-darwin13.4.0)
  CPU: Intel(R) Core(TM) i7-5557U CPU @ 3.10GHz
  WORD_SIZE: 64
  BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
  LAPACK: libopenblas64_
  LIBM: libopenlibm
  LLVM: libLLVM-3.3

Also works on:
Julia Version 0.5.1-pre+2
Commit f0d40ec (2016-09-20 03:34 UTC)
Platform Info:
  System: Darwin (x86_64-apple-darwin15.6.0)
  CPU: Intel(R) Core(TM) i7-5557U CPU @ 3.10GHz
  WORD_SIZE: 64
  BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
  LAPACK: libopenblas64_
  LIBM: libopenlibm
  LLVM: libLLVM-3.7.1 (ORCJIT, broadwell)

=#



# Store the hospitals in data structures to reduce the complexity of writing out all of the values
# to matrices and keeping track of all of the indices, etc.
#addprocs(2)
#include("/Users/austinbean/Desktop/dynhosp/Reboot.jl")


        ##### NB: Supply-side Data Creation Functions ######




"""
`MakeIt(Tex::EntireState, fip::Vector)`

Testing this:
MakeIt(ProjectModule.fips);
"""
function MakeIt(fip::Vector{Int64})
  Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
  for f in fip
    if f != 0
      push!(Tex.ms, Market( Array{hospital,1}(), Dict{Int64, hospital}(), f, Dict{Int64, Bool}()))
    end
  end
  Tex.mkts = Dict(m.fipscode => m for m in Tex.ms)
  return Tex
end



"""
`TXSetup(Tex::EntireState, data::Matrix; ...)`
 Takes an entire state and adds data from the imported choices returns a record with
 fipscodes containing hospitals with mostly empty field values

 Testing:
 Tex = MakeIt(ProjectModule.fips);
 TXSetup(Tex, ProjectModule.alldists);
OR:
 Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists);

"""
function TXSetup(Tex::EntireState, data::Matrix;
                 fidcol::Int64 = 4,
                 latcol::Int64 = 21,
                 longcol::Int64 = 22,
                 namecol::Int64 = 5,
                 fipscol::Int64 = 7,
                 intensivecol::Int64 = 11,
                 intermediatecol::Int64 = 20)
  for i = 1:size(data,1)
    fips = data[i, fipscol]
    if fips != 0
      level = 0
      if (data[i, intensivecol] == 1)&(data[i,intermediatecol]==0)
        level = 3
      elseif (data[i, intensivecol] == 0)&(data[i,intermediatecol]==1)
        level = 2
      else
        level = 1
      end
      availfids = FindFids(Tex.mkts[fips])
      if !in(data[i, fidcol], availfids)
        push!(Tex.mkts[fips].config,
        hospital( data[i, fidcol], #fid
                  data[i,latcol], #lat
                  data[i, longcol], #long
                  data[i, namecol], #name
                  fips,
                  level, # level
                  Array{Int64,1}(),
                  DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                  DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                  WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                  WeightVec([0.0]), #choice probs
                  Array{Float64,1}(),
                  neighbors(0,0,0,0,0,0,0,0,0), #neighbors
                  Array{Int64,1}(), # hood (array of fids)
                    0    , # beds added later.
                  false ) ) # perturbed or not?
      end
    end
    # push all hospital fid/ fips pairs into the directory.
    Tex.fipsdirectory[data[i, fidcol]] = fips # now for the whole state I can immediately figure out which market a hospital is in.
  end
  NeighborFix(Tex)
  InitChoice(Tex)
  ExpandDict(Tex)
  return Tex
end

"""
`FindFids(m::Market)`
Return all the fids in the market.
This is necessary in TXSetup to make sure that we only add each hospital to the market once.

Testing:
  Setup -
Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists);
  Testing:
FindFids(Tex.mkts[48453])
FindFids(Tex.mkts[48201])
"""
function FindFids(m::Market)
  outp::Array{Int64,1} = Array{Int64,1}()
  for el in m.config
    push!(outp, el.fid)
  end
  return outp
end

"""
`InitChoice(Tex::EntireState)`
Takes a newly created state and fixes all of the choice probabilities.

Testing:
Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists);
InitChoice(Tex)

#NB - InitChoice is already called in Tex.  Test with:
Tex.mkts[48453].config[1].chprobability
"""
function InitChoice(Tex::EntireState)
  for mk in keys(Tex.mkts)
    for el in Tex.mkts[mk].config
  #    levl=(-1,-1)
      if el.level == 1
        levl = (0,0)
        levels = MktSize(el.neigh)
        el.chprobability = WeightVec(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
      elseif el.level == 2
        levl = (1,0)
        levels = MktSize(el.neigh)
        el.chprobability = WeightVec(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
      elseif el.level == 3
        levl = (0,1)
        levels = MktSize(el.neigh)
        el.chprobability = WeightVec(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
      end
#      levels = MktSize(el.neigh)
#      el.chprobability = WeightVec(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
    end
  end
end




"""
`ExpandDict(Tex::EntireState)`
 Expand the market dictionaries so that they are filled with the hospitals

 NB: Can I run this within TXSetup??
"""
function ExpandDict(Tex::EntireState)
  for el in Tex.ms
    el.collection = Dict(i.fid => i for i in el.config)
    el.noneqrecord = Dict(i.fid => false for i in el.config)
    # I would like to append to each hospital a list of the others in the market, if it isn't already there.
    for hosp in el.config
      for hosp2 in el.config
        if hosp.fid != hosp2.fid
          if !in(hosp2.fid, hosp.hood) #TODO - is this condition necessary?  When will this be tripped?
            push!(hosp.hood, hosp2.fid)
          end
        end
      end
    end
  end
end



"""
`MakeNew(fi::Vector, dat::Matrix)`
Call this and the whole state with all markets should be created.
Should be called on "fips" or ProjectModule.fips and ProjectModule.alldists.
MakeNew(ProjectModule.fips, ProjectModule.alldists)
"""
function MakeNew(fi::Vector, dat::Matrix)
  return TXSetup(MakeIt(fi), dat)
end



# TODO: How are these functions different?
"""
`CreateEmpty(fi::Vector, dat::Matrix)`
 This creates an empty entire state record for the perturbed simulation.
"""
function CreateEmpty(fi::Vector, dat::Matrix)
  Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
  MakeIt(Tex, fi)
  TXSetup(Tex, dat)
  ExpandDict(Tex)
  return Tex
end

#NB: Creates hospital datastructure
#Texas = MakeNew(fips, ProjectModule.alldists );

      #### NB:  Supply-side  Printing Utilities to Display Simulation Outcomes




"""
`MarketPrint(mkt::Market)`
Prints the elements of the market record: name, neighbors, choice probabilities.
For debugging purposes to make sure things look right.
"""
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




"""
`NeighborsPrint(mkt::Market)`
Prints the name and neighbors of every hospital in the market `mkt`
"""
function NeighborsPrint(mkt::Market)
  for el in mkt.config
    println(el.name, " ", el.neigh)
  end
end




"""
`FacPrint(hosp::hospital)`
Prints out a hospital facility record: name, fid, fips, level, neighbors.
"""
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




"""
`NewEntrantLocation(mkt::Market)`
Takes the market, takes the mean location of all hospitals, adds normal noise to it.  ≈ 6 miles perturbation from mean.
"""
function NewEntrantLocation(mkt::Market)
  meanlat::Float64 = 0.0
  meanlong::Float64 = 0.0
  for el in mkt.config # over hospitals
    meanlat += el.lat
    meanlong += el.long
  end
  return [meanlat/size(mkt.config, 1) + rand(Normal(0, 0.1), 1)[1], meanlong/size(mkt.config, 1) + rand(Normal(0, 0.1), 1)[1]]::Array{Float64,1}
end



"""
`MktSize(n::Neighbors)`
Operates on a type `n` set of neighbors - returns the number of hospitals at levels 1, 2 and 3 across the distance categories
"""
function MktSize(n::neighbors)
  # takes a set of neighbors and returns the sum of levels 1, 2, 3 at the various distances.
  sum1::Int64 = n.level105 + n.level1515 + n.level11525
  sum2::Int64 = n.level205 + n.level2515 + n.level21525
  sum3::Int64 = n.level305 + n.level3515 + n.level31525
  return sum1, sum2, sum3
end


"""
`ChoicesAvailable(h::hospital)`
Takes a hospital, returns the choices available at that level as vectors.
e.g., returns the numbered choices available to it.
"""
function ChoicesAvailable{T<:Fac}(h::T)
  if h.level == 1
    return [10 2 1 11]::Array{Int64,2}
  elseif h.level == 2
    return [5 10 6 11]::Array{Int64,2}
  elseif h.level == 3
    return [4 3 10 11]::Array{Int64,2}
  else # exited
    return [-999 -999 -999 -999]::Array{Int64,2}
  end
end




"""
`LevelFunction(h::hospital, choice::Int64)`
Takes a hospital record and a choice and returns the corresponding level.
That is, what will the level be next period.
"""
function LevelFunction{T<:Fac}(h::T, choice::Int64)
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




"""
`NeighborAppend{T<:Fac}(elm::T, entrant::T)`
Takes two hospital records, computes the distance between them and adds a 1 to the relevant record in the neighborhood type.
Appends it to the hood of elm, which is a list of fids.  So this adds to both elm.neigh and elm.hood.
It is not symmetric - it appends entrant to elm, not vice versa.  Extended to include arguments of `chospital` type for
the counterfactual simulation.
"""
function NeighborAppend{T<:Fac}(elm::T, entrant::T)
  dist::Float64 = distance(elm.lat, elm.long, entrant.lat, entrant.long )
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




"""
`NeighborRemove{T<:Fac}(elm::T, entrant::T)`
takes two hospital records, computes the distance between them and subtracts 1 from the relevant record in the neighborhood type.
It removes the record of entrant FROM the record of elm.  Also not symmetric - removes entrant from elm's records, not the reverse.
"""
function NeighborRemove{T<:Fac}(elm::T, entrant::T)
  dist::Float64 = distance(elm.lat, elm.long, entrant.lat, entrant.long )
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




"""
`NeighborClean(state::EntireState)`
Takes every hospital in every market and sets all neighbors to 0
Also sets the distance/category "neighbors" to a vector of zeros.
"""
function NeighborClean(state::EntireState)
  for mkt in state.ms
    for hosp in mkt.config
      hosp.neigh::neighbors = neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0)
      hosp.hood::Array{Int64,1} = Array{Int64,1}()
    end
  end
end




"""
`NeighborFix(state::EntireState)`
For every hospital in the state, append all other hospitals within 25 miles, ignoring county boundaries.
Compare to StrictCountyNeighborFix which respects county boundaries.

"""
function NeighborFix(state::EntireState)
  for mk1 in state.ms
    for mk2 in state.ms
      for h1 in mk1.config
        for h2 in mk2.config
          if (h1.fid != h2.fid)&&(distance(h1.lat, h1.long, h2.lat, h2.long) < 25.0)
            NeighborAppend(h1, h2)
            NeighborAppend(h2, h1)
          end
        end
      end
    end
  end
end
 #=
 #variant
function NeighborFix(state::EntireState)
  for mkt1 in state.ms
    for h1 in mkt1.config
      for h2 in mkt1.config
        if h1.fid != h2.fid
          NeighborAppend(h1, h2)
        end
      end
    end
    # Now do all other counties in the state.
    for mkt2 in state.ms
      if mkt1.fipscode != mkt2.fipscode
        for hos1 in mkt1.config
          for hos2 in mkt2.config
            if distance(hos1.lat, hos1.long, hos2.lat, hos2.long) < 25
              NeighborAppend(hos1, hos2)
              NeighborAppend(hos2, hos1) # The function is not symmetric, so the second call is important.
            end
          end
        end
      end
    end
  end
end
  =#





"""
`StrictCountyNeighborFix(state::EntireState)`
 For every hospital in the state, append all other hospitals within 25 miles AND in the same county.
 More restrictive than NeighborFix.
"""
function StrictCountyNeighborFix(state::EntireState)
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




"""
`HospFindFirst(mkt::Market, hosp::hospital)`
looks for a hospital given by hosp in the market mkt, by searching for the fid and returning the index.
This is to search for the index of the hosp record in the mkt.config array.
"""
function HospFindFirst(mkt::Market, hosp::hospital)
  found = 0
  for el in 1:size(mkt.config,1)
    if mkt.config[el].fid == hosp.fid
      found = el
    end
  end
  return found
end




"""
`FidFindFirst(mkt::Market, fid::Int64)`
looks for a fid in the market, then returns the index of the fid in the mkt.config array
"""
function FidFindFirst(mkt::Market, fid::Int64)
  found = 0
  for el in 1:size(mkt.config,1)
    if mkt.config[el].fid == fid
      found = el
    end
  end
  return found
end




"""
`MarketCleaner(mkt::Market)`
Takes a whole market and then removes the records of any entrants.
"""
function MarketCleaner(mkt::Market)
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




"""
`HospUpdate(hosp::hospital, choice::Int; update = false)`
Takes a hospital record and updates the probabilities of the choices.
"""
function HospUpdate{T<:Fac}(hosp::T, choice::Int; update = false)
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




"""
`HospPerturb(hosp::hospital, choice::Int, eps::Float64)`
Takes a hospital record and updates the probabilities of the choices.
and then perturbs them using the perturb function.
"""
function HospPerturb(hosp::hospital, choice::Int, eps::Float64)
  # TODO: Is there a change which needs to be made when other firms do something weird?  Does that make sense?
  levl = (-1, -1)
 if (hosp.level!=choice) # want to be able to force this to rerun when the data is cleaned again.
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
     prs = perturb(prs, eps, false)
     return WeightVec(vec(prs))
   else # choice = -999
     return WeightVec([1.0]) #TODO: one option, no choices ??  Might need four options [1.0 1.0 1.0 1.0]
   end
 else
   return hosp.chprobability
  end
end

      ##### NB: Demand-side Data Structure Creation.


"""
`function CreateZips(alld::Array,Tex::EntireState)`
Now the correct distances are being imported as ProjectModule.alldists.
But in this version do check to make sure that the distance between any fac and zip is
less than 25 miles.
Note that this indexes into the Array according to the column zipcol.
Create the state first with:
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, XXXXX);
CreateZips(ProjectModule.alldists, Tex)
"""
function CreateZips(alld::Array,
                     Tex::EntireState;
                     zipcol::Int64 = 1,
                     fidcol::Int64 = 4,
                     bedcol::Int64 = 10,
                     latcol::Int64 = 2,
                     longcol::Int64 = 3,
                     dat::Array{Any,2} = ProjectModule.alldists,
                     datfidloc::Int64 = 4,
                     bedmean::Float64 = round(mean(alld[:,bedcol])))
  ppatients::patientcollection = patientcollection( Dict{Int64, zip}() )
  # unfound = Array{Int64,1}()
  # found = Array{Int64,1}()
  # names = Array{AbstractString, 1}()
  fids = unique(dat[:,fidcol])
  matched::Int64 = 0
  ppatients.zips = Dict(k=> zip(k, 0, Dict{Int64,ProjectModule.Fac}(), Dict{Int64,Float64}(),                                                                              # zipcode, public health region, facilities, hospital FE's.
                           Dict{Int64,Float64}(), Dict{Int64, Float64}(),                                                                                     # private det utilities, medicaid det utilities.
                           0.0, 0.0,
                           coefficients(ProjectModule.privatedistance_c, ProjectModule.privatedistsq_c, ProjectModule.privateneoint_c, ProjectModule.privatesoloint_c, ProjectModule.privatedistbed_c, ProjectModule.privateclosest_c),
                           coefficients(ProjectModule.medicaiddistance_c, ProjectModule.medicaiddistsq_c, ProjectModule.medicaidneoint_c, ProjectModule.medicaidsoloint_c, ProjectModule.medicaiddistbed_c, ProjectModule.medicaidclosest_c),     # lat, long, private coefficients, medicaid coefficients
                           patientcount(0,0,0,0,0,0,0), patientcount(0,0,0,0,0,0,0)) for k in unique(alld[:,zipcol]))
  for row = 1:size(alld,1)
    if alld[row,fidcol] != 0
      if in(alld[row, fidcol], fids)
        #TODO - ADD THE LAT AND LONG CODE HERE.
        ppatients.zips[alld[row,zipcol]].facilities[alld[row,fidcol]] = Tex.mkts[ Tex.fipsdirectory[alld[row, fidcol]] ].collection[alld[row,fidcol]]
        if alld[row,bedcol] != 0
          Tex.mkts[Tex.fipsdirectory[alld[row, fidcol]]].collection[ alld[row,fidcol]].bedcount = alld[row, bedcol]
        else
           Tex.mkts[Tex.fipsdirectory[alld[row, fidcol]]].collection[ alld[row,fidcol]].bedcount = bedmean
        end
        if ppatients.zips[alld[row,zipcol]].lat == 0
          ppatients.zips[alld[row,zipcol]].lat = alld[row, latcol]
          ppatients.zips[alld[row,zipcol]].long = alld[row, longcol]
        end
        # push!(found, alld[row, fidcol])
      else
        # push!(unfound, alld[row, fidcol])
        # push!(names, alld[row, :facility])
      end
    end
  end
  return  ppatients   #, unique(unfound), unique(names), unique(found)
end

"""
`AddOO(patients::patientcollection)`
Add an outside option with zero utility to the choices of patients.
"""
function AddOO(pats::patientcollection)
  for el in keys(pats.zips)
    pats.zips[el].pdetutils[0] = 0.0
    pats.zips[el].mdetutils[0] = 0.0
  end
end



"""
`FillPPatients(pats::patientcollection, imported::Matrix; ziploc = 101, drgloc = 104)`
Takes the imported matrix of *privately-insured* patients and records the number at each DRG 385-391 in each zip record.
There is a separate function for the Medicaid patients.
There is now one zip which patients are in but which is not in pats as constructed above.
"""
function FillPPatients(pats::patientcollection, imported::Matrix;
                       ziploc::Int64 = 101,
                       drgloc::Int64 = 104,
                       notavail::Array{Float32,1} = setdiff(imported[:,ziploc] ,keys(pats.zips)))
  notfound = Array{Int64,1}()
  for row in 1:size(imported, 1)
    if !in(imported[row,ziploc], notavail)
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
  end
  return pats
end




"""
`FillMPatients(pats::patientcollection, imported::Matrix; ziploc = 101, drgloc = 104)`
Takes the imported matrix of *Medicaid* patients and records the number at each DRG 385-391 in each zip record.
There is a separate function for the privately-insured patients.
"""
function FillMPatients(pats::patientcollection, imported::Matrix;
                       ziploc::Int64 = 101,
                       drgloc::Int64 = 104,
                       notavail::Array{Float32,1} = setdiff(imported[:,ziploc] ,keys(pats.zips)))
  notfound = Array{Int64,1}()
  for row in 1:size(imported, 1)
    if !in(imported[row,ziploc], notavail)
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
  end
  return pats;
end




"""
`FillPPatients(pats::patientcollection, private::Matrix, medicaid::Matrix)`
 Adds the privately insured and medicaid patients to the zip records.
"""
function FillPatients(pats::patientcollection, private::Matrix, medicaid::Matrix)
  pats = FillPPatients(pats, private)
  pats = FillMPatients(pats, medicaid)
end


"""
`CheckPats(pats::patientcollection)`
This will compute the total sum of all patients in all zips in the patientcollection.
This is just to check that all are being added as expected.
"""
function CheckPats(pats::patientcollection)
  count::Int64 = 0
  for k in keys(pats.zips)
    count += sum(pats.zips[k].mpatients)
    count += sum(pats.zips[k].ppatients)
  end
  return count
end



"""
`ComputeDetUtil(zipc::zip, fid::Int64, p_or_m::Bool)`
Computes the deterministic component of utility for each hospital in the zip "zipc".
Maps exited facilites to have deterministic utility -999
Works on private and medicaid patients by setting p_or_m to true or false, respectively.
Has been written to accomodate hospital FE's when available.
"""
function ComputeDetUtil(zipc::zip, fid::Int64, p_or_m::Bool)
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





"""
`UpdateDeterministic(collect::patientcollection)`
Computes the deterministic component of the utility - updates every firm every time it is called.
Is called during the Eq and Non-eq simulations.
"""
function UpdateDeterministic(collect::patientcollection)
  for el in keys(collect.zips) #iterates over zips
    for fid in keys(collect.zips[el].facilities) # iterates over dict of facilities within zip.
      collect.zips[el].mdetutils[fid] = ComputeDetUtil(collect.zips[el], fid, false)
      collect.zips[el].pdetutils[fid] = ComputeDetUtil(collect.zips[el], fid, true)
    end
  end
end




"""
`NewPatients(Tex::EntireState; dists = ProjectModule.alldists, phrloc = 103, pins = pinsured, pmed = pmedicaid)`
this creates the whole collection of patients.  0.7 seconds.  Pretty slow.
It must take an existing EntireState record to link the hospitals.

#TODO - there are patients who don't get added.  Their zips are not in the list but are in the inpatient discharge.
Texas = MakeNew(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);
"""
function NewPatients(Tex::EntireState;
                     dists::Array{Any,2} = ProjectModule.alldists,
                     phrloc::Int64 = 103,
                     pins = ProjectModule.pinsured,
                     pmed = ProjectModule.pmedicaid)
  patients = CreateZips(dists, Tex) #NB: This needs to take the whole state so that the hosps in zips point to the same underlying record.
  patients = FillPatients(patients, pins, pmed)
  UpdateDeterministic(patients)
  AddOO(patients)
  return patients
end

"""
`NewPatientsTest(;pins = ProjectModule.pinsured, pmed = ProjectModule.pmedicaid)`
Test to make sure that non-empty zips are created for each of the zips in the patient data.
What it does:
- For patients in the *data* - count the number at each zip code.
- For patients generated by NewPatients(), IF the number of patients at that zip is zero, but the number in the data isn't, prints a warning.
Call it:
Texas = MakeNew(ProjectModule.fips, ProjectModule.alldists);
NewPatientsTest(Texas)
"""
function NewPatientsTest(Tex::EntireState;
                     dists::Array{Any,2} = ProjectModule.alldists,
                     phrloc::Int64 = 103,
                     ziploc::Int64 = 101,
                     pins::Array{Float32,2} = ProjectModule.pinsured,
                     pmed::Array{Float32,2} = ProjectModule.pmedicaid)
  patients = CreateZips(dists, Tex) # creates empty zip codes.
  patients = FillPatients(patients, pins, pmed)
  UpdateDeterministic(patients)
  AddOO(patients)
  pzips::Dict{Int64, Int64} = Dict(k => 0 for k in unique(pins[:,ziploc]))
  for n in 1:size(pins, 1)
    if haskey(pzips, pins[n,ziploc])
      pzips[pins[n,ziploc]] += 1
    else
      println("No key for ", pins[n,ziploc])
    end
  end
  mzips::Dict{Int64, Int64} = Dict(k => 0 for k in unique(pmed[:,ziploc]))
  for n in 1:size(pmed,1)
    if haskey(mzips, pmed[n,ziploc])
      mzips[pmed[n,ziploc]] += 1
    else
      println("No key for ", pmed[n, ziploc])
    end
  end
  nothere::Array{Int64,1} = Array{Int64,1}()
  for el in keys(patients.zips)
    if (sum(patients.zips[el].ppatients) == 0)
      if haskey(pzips, el)
        if pzips[el] != 0
          println("At zip ", el)
          println("Private patients number ", sum(patients.zips[el].ppatients))
          println("But in the data: ", pzips[el])
        end
      else
        println("absent key ", el )
        push!(nothere, el)
      end
    end
    if (sum(patients.zips[el].mpatients) == 0)
      if haskey(mzips, el)
        if mzips[el] != 0
          println("At zip ", el)
          println("Medicaid patients number ", sum(patients.zips[el].mpatients))
          println("But in the data: ", mzips[el])
        end
      else
        println("absent key ", el)
        push!(nothere, el)
      end
    end
  end
  return unique(nothere)
end


"""
`FillState(Tex::EntireState)`
fills the entire state record with elements of the chospital type
for the counterfactual only.
Note that this needs to be called AFTER the function `CMakeIt(Tex::EntireState, fip::Vector)` is called
on an empty state record.
This version also adds all of the `chospitals` directly to the `Market.collection` dictionary.
Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists);
"""
function FillState(Tex::EntireState, data::Matrix;
                   fidcol::Int64 = 4,
                   latcol::Int64 = 21,
                   longcol::Int64 = 22,
                   namecol::Int64 = 5,
                   fipscol::Int64 = 7,
                   intensivecol::Int64 = 11,
                   intermediatecol::Int64 = 20)
  for i = 1:size(data,1)
    fips = data[i, fipscol] # fipscode
    if fips != 0
      level::Int64 = 0
      if (data[i, intensivecol] == 1)&(data[i,intermediatecol]==0)
        level = 3
      elseif (data[i, intensivecol] == 0)&(data[i,intermediatecol]==1)
        level = 2
      else
        level = 1
      end
      # Create the new record.
      availfids = FindFids(Tex.mkts[fips])
      if !in(data[i, fidcol], availfids)
        newh = ProjectModule.chospital( data[i, fidcol],
                  data[i,latcol],
                  data[i, longcol],
                  data[i, namecol],
                  fips,
                  level,
                  level,
                  Array{Int64,1}(), #volume
                  Array{Int64, 1}(), #mortality
                  Array{Float64,1}(), #ppayoff
                  Array{Float64,1}(), #mpayoff
                    0    , # beds added later.
                  LBW(0,0,0,0,0,0), # LBW Infants.
                  false, # has intensive
                  false )
        push!(Tex.mkts[fips].config,newh)
        Tex.mkts[fips].collection[data[i,fidcol]] = newh # finished.
      end
    end
    # push all hospital fid/ fips pairs into the directory.
    Tex.fipsdirectory[data[i, fidcol]] = fips # now for the whole state I can immediately figure out which market a hospital is in.
  end
  #TODO - add the zero hospital at the zero market.
  return Tex
end



"""
`CMakeIt(Tex::EntireState, fip::Vector)`
Perhaps poor practice to use Eval in this way, but generates markets named m*fipscode* for any fipscode in the vector fip.
"""
function CMakeIt(Tex::EntireState, fip::Vector)
  for el in fip
    if el != 0
      el = eval(parse("m$el = Market( Array{chospital,1}(), Dict{Int64, chospital}(), $el, Dict{Int64, Bool}())"))
      push!(Tex.ms, el)
    end
  end
  Tex.mkts = Dict(m.fipscode => m for m in Tex.ms)
  # Tex.mkts = [ m.fipscode => m for m in Tex.ms] # this is the pre0.5 generator syntax
end


"""
`SetLevel(mkt::Market, fid::Int64, level::Int64; special::Int64 = 3)`
This function should set all of the facility levels in a given market to `level`, except for that specified by fid.
That one is set to the optional argument `special`, which by default is 3.
"""
function SetLevel(mkt::Market, sfid::Int64, level::Int64; special::Int64 = 3)
  for el in mkt.config
    if el.fid != sfid
      el.level = level
    elseif el.fid == sfid
      el.level = special
    end
  end
end



"""
`FindUndone(mkt::Market)`
Takes a market, returns a vector of the fids which do not have "finished" set to true.
"""
function FindUndone(mkt::Market)
  outp = Array{Int64,1}()
  for el in mkt.config
    if !el.finished
        push!(outp, el.fid)
    end
  end
  return outp
end



    ### NB: Zip code record printing utility.




"""
`PrintZip(zi::zip)`
Prints the fid and the name of the facilities attached to the zips.
"""
function PrintZip(zi::zip)
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




"""
`WhichZips(pats::patientcollection, fid::Int64)`
Takes a patientcollection and tells me which zips have the hospital fid
"""
function WhichZips(pats::patientcollection, fid::Int64)
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



"""
`CalcWTP(zipc::zip)`
Takes the deterministic component of utility for the privately insured patients and returns a WTP measure.
Output is sent to WTPMap.
#NB - adding an outside option here.
"""
function CalcWTP(zipc::zip)
  outp = Dict(j=> 0.0 for j in keys(zipc.pdetutils))
#  outp[0] = 0.0 # maps 0, the OO fid, to 0.0, the OO utility.
  interim = 0.0
  for el in keys(zipc.pdetutils)
    interim +=  (outp[el] = exp(zipc.pdetutils[el]) )
  end
  return Dict( j=> outp[j]/interim for j in keys(outp))
end



"""
`WTPMap(pats::patientcollection, Tex::EntireState)`
Takes a patient collection and an entire state and returns a dict{fid, WTP}
computed by calling CalcWTP.  Right now it ignores Inf and NaN.
Input is from CalcWTP.  Output is sent to WriteWTP
# NB: changing the try-catch to haskey.
"""
function WTPMap(pats::patientcollection, Tex::EntireState)
  outp = Dict(j=>0.0 for j in keys(Tex.fipsdirectory))
  for zipc in keys(pats.zips)
    vals = CalcWTP(pats.zips[zipc])
    for el in keys(vals) # What to do about key errors?  there will be some.
      if haskey(outp, el)
        if vals[el] != 1 & !isnan(vals[el]) # there should be no way for this to be 1 anyway.
          outp[el] += log(1/(1-vals[el]))
        elseif (vals[el] == 1)||(isnan(vals[el]))
            println(zipc, "  ", vals[el])
        end
      else # key absent
        # facility missing.
      end
    end
  end
  return outp # gives a dict{fid, WTP} back
end




"""
`WriteWTP(reslt::Dict{Int64, Float64}, Tex::EntireState)`
Takes a dict of {fid, WTP} and writes it out by DRG.
Works on the output of WTPMap
"""
function WriteWTP(reslt::Dict{Int64, Float64}, Tex::EntireState)
  for els in keys(reslt)
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w385, reslt[els])
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w386, reslt[els])
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w387, reslt[els])
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w388, reslt[els])
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w389, reslt[els])
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w390, reslt[els])
    push!(Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w391, reslt[els])
  end
end



"""
`GenPChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))`
The patient choice is max \bar{U} + ϵ, but we have \bar{U} from Compute Det Util and we know how many patients there are in
the privately insured category from FillPPatients.  This returns a dict of fids and patient counts, where patient counts are
generated by repeatedly finding max i = 1, ..., N \bar{U}_i + ϵ_i.  Note the corresponding GenMChoices below.
"""
function GenPChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  outp = Dict( j=> patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory) )
  # outp = [ j => patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory)] # output is a {FID, patientcount} dictionary. This is pre 0.5 syntax.
  for zipcode in keys(pats.zips)
    if pats.zips[zipcode].pdetutils.count > 0
      utils = hcat([ [k1,pats.zips[zipcode].pdetutils[k1]] for k1 in keys(pats.zips[zipcode].pdetutils)]...)
      temparr = zeros(size(utils, 2))
      for k = 1:pats.zips[zipcode].ppatients.count385
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count385 += 1
      end
      for k = 1:pats.zips[zipcode].ppatients.count386
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count386 += 1
      end
      for k = 1:pats.zips[zipcode].ppatients.count387
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count387 += 1
      end
      for i=1:pats.zips[zipcode].ppatients.count388
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count388 += 1
      end
      for i = 1:pats.zips[zipcode].ppatients.count389
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count389 += 1
      end
      for i=1:pats.zips[zipcode].ppatients.count390
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count390 += 1
      end
      for i = 1:pats.zips[zipcode].ppatients.count391
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count391 += 1
      end
    end
  end
  return outp
end




"""
`GenMChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))`
The patient choice is max \bar{U} + ϵ, but we have \bar{U} from Compute Det Util and we know how many patients there are in
the Medicaid category from FillMPatients.  This returns a dict of fids and patient counts, where patient counts are
generated by repeatedly finding max i = 1, ..., N \bar{U}_i + ϵ_i.  Note the corresponding GenPChoices above.
"""
function GenMChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
   outp = Dict( j=> patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory) )
  #  outp = [ j => patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory)] # output is a {FID, patientcount} dictionary.  This is the pre0.5 syntax
  for zipcode in keys(pats.zips)
    if pats.zips[zipcode].mdetutils.count > 0
      utils = hcat([ [k1,pats.zips[zipcode].mdetutils[k1]] for k1 in keys(pats.zips[zipcode].mdetutils)]...)
      temparr = zeros(size(utils, 2))
      for k = 1:pats.zips[zipcode].mpatients.count385
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count385 += 1
      end
      for k = 1:pats.zips[zipcode].mpatients.count386
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count386 += 1
      end
      for k = 1:pats.zips[zipcode].mpatients.count387
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count387 += 1
      end
      for i=1:pats.zips[zipcode].mpatients.count388
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count388 += 1
      end
      for i = 1:pats.zips[zipcode].mpatients.count389
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count389 += 1
      end
      for i=1:pats.zips[zipcode].mpatients.count390
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count390 += 1
      end
      for i = 1:pats.zips[zipcode].mpatients.count391
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count391 += 1
      end
    end
  end
  return outp
end




"""
`PHistoryAdd(hos::hospital, cnt::patientcount)`
Maps patientcount to the private demand history
"""
function PHistoryAdd(hos::hospital, cnt::patientcount)
  push!(hos.pdemandhist.demand385, cnt.count385)
  push!(hos.pdemandhist.demand386, cnt.count386)
  push!(hos.pdemandhist.demand387, cnt.count387)
  push!(hos.pdemandhist.demand388, cnt.count388)
  push!(hos.pdemandhist.demand389, cnt.count389)
  push!(hos.pdemandhist.demand390, cnt.count390)
  push!(hos.pdemandhist.demand391, cnt.count391)
end




"""
`MHistoryAdd(hos::hospital, cnt::patientcount)`
Maps patientcount to the Medicaid demand history.
"""
function MHistoryAdd(hos::hospital, cnt::patientcount)
  push!(hos.mdemandhist.demand385, cnt.count385)
  push!(hos.mdemandhist.demand386, cnt.count386)
  push!(hos.mdemandhist.demand387, cnt.count387)
  push!(hos.mdemandhist.demand388, cnt.count388)
  push!(hos.mdemandhist.demand389, cnt.count389)
  push!(hos.mdemandhist.demand390, cnt.count390)
  push!(hos.mdemandhist.demand391, cnt.count391)
end




"""
`PDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)`
Maps Private patient demand out to the state record.
"""
function PDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)
  for el in keys(patd)
    PHistoryAdd(Tex.mkts[Tex.fipsdirectory[el]].collection[el], patd[el])
  end
end




"""
`MDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)`
Maps Medicaid Patient demand out to the state record.
"""
function MDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)
  for el in keys(patd)
    MHistoryAdd(Tex.mkts[Tex.fipsdirectory[el]].collection[el], patd[el])
  end
end


"""
`PatientDraw(ppat::Dict, mpat::Dict; bins = collect(1:13), weightprobs = bwprobs, admitprobs = nicuprobs)`
Takes two dictionaries - one of medicaid and the other of private patients, adds all the patients together,
then for each patient draws a weight.  Those LBW are assumed to be admitted.  At higher weights admission
probs are drawn from the empirical distribution from NCHS birth certificate data.  The return type is a
`Dict{Int64, LBW}` dictionary of fid's and low birth weight types.
This assumes all LBW patients get admitted and then draws admit probs for those in higher categories.  Note
that not all are born in hospitals with NICU's, so this is a bit problematic.
This function also returns the two arguments ppat and mpat.  These are needed later.
"""

function PatientDraw(ppat::Dict, mpat::Dict, Tex::EntireState;
                     bins = collect(1:13),
                     weightpr = weightprobs[:,2],
                     admitprobs = nicuprobs[:,2],
                     w1 = WeightVec([1-admitprobs[1],admitprobs[1]]), # next few lines unused.
                     w2 = WeightVec([1-admitprobs[2],admitprobs[2]]),
                     w3 = WeightVec([1-admitprobs[3],admitprobs[3]]),
                     w4 = WeightVec([1-admitprobs[4],admitprobs[4]]),
                     w5 = WeightVec([1-admitprobs[5],admitprobs[5]]),
                     w6 = WeightVec([1-admitprobs[6],admitprobs[6]]),
                     w7 = WeightVec([1-admitprobs[7],admitprobs[7]]), # here and below all used.
                     w8 = WeightVec([1-admitprobs[8],admitprobs[8]]),
                     w9 = WeightVec([1-admitprobs[9],admitprobs[9]]),
                     w10 = WeightVec([1-admitprobs[10],admitprobs[10]]),
                     w11 = WeightVec([1-admitprobs[11],admitprobs[11]]),
                     w12 = WeightVec([1-admitprobs[12],admitprobs[12]]),
                     w13 = WeightVec([1-admitprobs[13],admitprobs[13]]))
  outp::Dict{Int64, LBW} = Dict( k=>LBW(0,0,0,0,0,0) for k in keys(Tex.fipsdirectory)) # empty dictionary of fids/LBW record types
  for el in keys(ppat)
    totl = sum(ppat[el] + mpat[el])
    for i = 1:totl
      bwt = sample(bins, WeightVec(weightpr)) # sample from the distribution of birthweights.
      if bwt == 1                                                               # all at this weight are being admitted
        outp[el].bt05 += 1
      elseif bwt == 2                                                           # all at this weight are being admitted
        outp[el].bt510 += 1
      elseif bwt == 3                                                           # all at this weight are being admitted
        outp[el].bt510 += 1
      elseif bwt == 4                                                           # all at this weight are being admitted
        outp[el].bt1015 += 1
      elseif bwt == 5                                                           # all at this weight are being admitted
        outp[el].bt1015 += 1
      elseif bwt == 6                                                           # all at this weight are being admitted
        outp[el].bt1520 += 1
      elseif bwt == 7                                                           # all at this weight are being admitted
        nicuadmit = sample([0,1], w7)
        if nicuadmit == 1
          outp[el].bt2025 += 1                                                  # at this weight admissions are stochastic
        end
      elseif bwt == 8
        nicuadmit = sample([0,1], w8)
        if nicuadmit == 1
          outp[el].bt2580 += 1                                                  # at this weight admissions are stochastic
        end
      elseif bwt == 9
        nicuadmit = sample([0,1], w9)
        if nicuadmit == 1
          outp[el].bt2580 += 1                                                  # at this weight admissions are stochastic
        end
      elseif bwt == 10
        nicuadmit = sample([0,1], w10)
        if nicuadmit == 1
          outp[el].bt2580 += 1                                                  # at this weight admissions are stochastic
        end
      elseif bwt == 11
        nicuadmit = sample([0,1], w11)
        if nicuadmit == 1
          outp[el].bt2580 += 1                                                  # at this weight admissions are stochastic
        end
      elseif bwt == 12
        nicuadmit = sample([0,1], w12)
        if nicuadmit == 1
          outp[el].bt2580 += 1                                                  # at this weight admissions are stochastic
        end
      else #bwt == 13
        nicuadmit = sample([0,1], w13)
        if nicuadmit == 1
          outp[el].bt2580 += 1                                                  # at this weight admissions are stochastic
        end
      end
    end
  end
  return outp, ppat, mpat
end




"""
`AllMortality(d::Dict{Int64, LBW})`
The return of `PatientDraw` is a dictionary of {fid, LBW}.  Take that volume and convert it to a mortality
rate.  Then apply the mortality rate to the LBW record.  The elements keys(d) will be fids.
"""
function AllMortality(d::Dict{Int64, LBW}, Tex::EntireState)
  outp = Dict(j => 0 for j in keys(d))
  for el in keys(d)
    outp[el] = floor(sum(d[el])*VolMortality(sum(d[el]), Tex.mkts[ Tex.fipsdirectory[el] ].collection[el].level))         # Function calls the level too
  end
  return outp
end



"""
`VolMortality(v::Int64)`
This function needs to return a mortality rate for the volume.  That is, take the number of patients and return the mortality
rate as a function of the patient volume.  This data comes from Chung, Phibbs, Boscardin, et al Medical Care 2010
"The Effect of Neonatal Intensive Care Level and Hospital Volume on Mortality of Very Low Birth Weigh Infants"
"""


function VolMortality{T<:Real}(v::T, lev::Int64)
  if lev == 1
    if (v>=0)&(v<10)
      return (-0.72*v + 19.9)/100
    elseif (v>=10)&(v<26)
       return (-0.72*v + 19.9)/100
    else # rate ceases to decline after 10.
      return 0.127
    end
  elseif lev == 2
    if (v>=0)&(v<10)
      return (-0.46*v + 19.16)/100
    elseif (v>=10)&(v<25)
      return (-0.053*v + 15.03)/100
    elseif (v>=25)&(v<283)
      return 0.13705 # Rate ceases to decline after 25.  alternatively: (-0.053*v+15.03)/100
    else
      return 0.13705
    end
  else lev == 3
    if (v>= 0)&(v<25)
      return (-0.052*v + 16.95)/100
    elseif (v>=25)&(v<50)
      return (-0.034*v + 16.5)/100
    elseif (v>=50)&(v<100)
      return (-0.046*v + 17.1)/100
    elseif (v>=100)&(v<371)
      return 0.125 # the rate ceases to decline after 100. alternatively: (-0.046*v + 17.1)/100
    else
      return 0.125  #
    end
  end
end


    ### NB: The business of the simulation.




"""
`NewSim(T::Int, Tex::EntireState, pats::patientcollection; entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002] )`
Runs a T period simulation using the whole state and whole collection of patient records.
"""
function NewSim(T::Int, Tex::EntireState, pats::patientcollection; entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002] )
  for i = 1:T
    WriteWTP(WTPMap(pats, Tex), Tex)
    PDemandMap(GenPChoices(pats, Tex), Tex)
    MDemandMap(GenMChoices(pats, Tex), Tex)
    for el in Tex.ms
      entrant = sample(entrants, WeightVec(entryprobs))
      if entrant != 0
        entloc = NewEntrantLocation(el)                                                        # called on the market
        newfid = -floor(rand()*1e6)-1000000                                                    # all entrant fids negative to facilitate their removal later.
        entr = hospital( newfid, entloc[1], entloc[2], " Entrant $newfid ", el.fipscode, entrant, [entrant],
                         DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                         DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                         WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                         WeightVec([0.1, 0.1, 0.1, 0.1]), Array{Float64,1}(), neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0), Array{Int64, 1}(), 0, false)
        push!(el.config, entr)                                                                 # need to create a new record for this hospital in the market
        el.collection[newfid] = entr
        for elm in el.config                                                                   # need to add it to the dictionary too:
          NeighborAppend(elm, entr)
          NeighborAppend(entr, elm)
        end
         HospUpdate(entr, entrant)                                                             #"entrant" is the level
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
    UpdateDeterministic(pats)                                                                  # Updates deterministic component of utility
  end
  return Tex                                                                                   # Returns the whole state so the results can be written out.
end

# Texas = MakeNew(fips, ProjectModule.alldists);
# patients = NewPatients(Texas);

#  Tex2 = NewSim(3, Texas, patients);
#  EmpTex = CreateEmpty(fips, XXX XXX );




"""
`Termination(EmTex::EntireState)`
Takes an entire state (or the empty state for data recording) and returns "true" when every facility has been perturbed.
"""
function Termination(EmTex::EntireState)
  isdone = true
  for mark in keys(EmTex.mkts) # iterates over markets
    isdone = (isdone)&(reduce(&, [ EmTex.mkts[mark].noneqrecord[i] for i in keys(EmTex.mkts[mark].noneqrecord) ] ))
  end
  return isdone
end




"""
`PSim(T::Int64 ; di = ProjectModule.alldists, fi = fips, entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002])`
Runs a perturbed simulation - for each market, while there are hospitals I have not perturbed, runs a sim with one perturbed and the rest not.
The results are stored in EmptyState, which is an EntireState record instance.
"""
function PSim(T::Int64 ; di = ProjectModule.alldists, fi = ProjectModule.fips, entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002])  # fi = fips,
  EmptyState = CreateEmpty(fi, di);                                                                     # This is just a container of EntireState type - does not need linking.
  termflag = true                                                                                       # Initializes the termination flag.
  counter = 1
  while termflag                                                                                        # true if there is some hospital which has not been perturbed.
    currentfac = Dict{Int64, Int64}()                                                                   # Dict{FID, fipscode} = {key, value}
    Tex = MakeNew(fi, di);                                                                              # NB: New state every time - this is kind of inefficient.
    pats = NewPatients(Tex);                                                                            # NB: New patient collection, linked to the new state.  Must be created AFTER "Tex."
    for el in keys(EmptyState.mkts)
      if !reduce(&, [ EmptyState.mkts[el].noneqrecord[i] for i in keys(EmptyState.mkts[el].noneqrecord)])
        pfids = prod(hcat( [ [i, !EmptyState.mkts[el].noneqrecord[i]] for i in keys(EmptyState.mkts[el].noneqrecord) ]...) , 1)
        pfid = pfids[findfirst(pfids)]                                                                  # takes the first non-zero element of the above and returns the element.
        currentfac[EmptyState.fipsdirectory[pfid]] = pfid                                               # Now the Key is the fipscode and the value is the fid.
        for hos in Tex.mkts[el].config
          if hos.fid == pfid
            hos.perturbed = true
          else
            hos.perturbed = false
          end
        end
      end
    end
    pmarkets = unique(keys(currentfac))                                                                # picks out the unique fipscodes remaining to be done.
    for i = 1:T
      WriteWTP(WTPMap(pats, Tex), Tex)
      PDemandMap(GenPChoices(pats, Tex), Tex)
      MDemandMap(GenMChoices(pats, Tex), Tex)
      for el in Tex.ms
        if in(el.fipscode, pmarkets) #NB: in( collection, element) !!
          entrant = sample(entrants, WeightVec(entryprobs))
          if entrant!= 0
            entloc = NewEntrantLocation(el)                                                            # called on the market
            newfid = -floor(rand()*1e6)-1000000                                                        # all entrant fids negative to facilitate their removal.
            entr = hospital( newfid, entloc[1], entloc[2], " Entrant $newfid ", el.fipscode, entrant, [entrant],
                             DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                             DemandHistory( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                             WTP( Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}(),  Array{Float64,1}(), Array{Float64,1}() ),
                             WeightVec([0.1, 0.1, 0.1, 0.1]), Array{Float64,1}(), neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0), Array{Int64, 1}(), 0, true) # entrants never perturbed.
                             push!(el.config, entr)                                                   # need to create a new record for this hospital in the market
             el.collection[newfid] = entr                                                             # need to add it to the dictionary too:
             for elm in el.config
               NeighborAppend(elm, entr)
               NeighborAppend(entr, elm)
             end
             HospUpdate(entr, entrant) #entrant is the level
          end
          for elm in el.config
             if !elm.perturbed                                                                        # not perturbed, i.e., "perturbed" == false
               action = sample( ChoicesAvailable(elm), elm.chprobability )                            # Take the action
               push!(elm.probhistory, elm.chprobability[ findin(ChoicesAvailable(elm), action)[1] ])  # Record the prob with which the action was taken.
               newchoice = LevelFunction(elm, action)                                                 # What is the new level?
               elm.chprobability = HospUpdate(elm, newchoice)                                         # What are the new probabilities, given the new level?
               elm.level = newchoice                                                                  # Set the level to be the new choice.
               push!(elm.levelhistory, newchoice)
             else # perturbed.
               action = sample( ChoicesAvailable(elm), HospPerturb(elm, elm.level,0.05))
               push!(elm.probhistory, elm.chprobability[findin(ChoicesAvailable(elm), action)[1]])
               newchoice = LevelFunction(elm, action)
               elm.chprobability = HospUpdate(elm, newchoice)
               elm.level = newchoice
               push!(elm.levelhistory, newchoice)
             end
          end
        end
      end
      UpdateDeterministic(pats)                                                                       # Updates deterministic component of utility for all patients and zips.
    end
    fipst = 0; fidt = 0;
    for fips in pmarkets                                                                              # a collection of fips codes
      EmptyState.mkts[fips].collection[ Tex.mkts[fips].collection[currentfac[fips]].fid ] = Tex.mkts[fips].collection[currentfac[fips]]
      EmptyState.mkts[fips].noneqrecord[ Tex.mkts[fips].collection[currentfac[fips]].fid ] = true     # update the value in the non-equilibrium record sim to true.
      for num in 1:size(EmptyState.mkts[fips].config,1)                                               # iterate over the market config, which is an array.
        if EmptyState.mkts[fips].config[num].fid == Tex.mkts[fips].collection[currentfac[fips]].fid   # check for equality in the fids
          EmptyState.mkts[fips].config[num] = Tex.mkts[fips].collection[currentfac[fips]]
        end
      end
    end
    termflag = !Termination(EmptyState)                                                               # Checks the termination condition over every market in the state.
    counter += 1
  end # of while
  return EmptyState
end


#  Perturbed = PSim(10);




"""
`TransitionGen(current::Int64, previous::Int64)`
Generates counts of transitions - checks current level against previous.
"""
function TransitionGen(current::Int64, previous::Int64)
  transitions = zeros(Int64, 9)
  if  (previous==1)&(current==2)
    transitions[1] += 1
  elseif (previous==1)&(current==3)
    transitions[2] += 1
  elseif (previous==1)&(current==-999)
    transitions[3] += 1
  elseif (previous==2)&(current==1)
    transitions[4] += 1
  elseif (previous==2)&(current==3)
    transitions[5] += 1
  elseif (previous==2)&(current==-999)
    transitions[6] += 1
  elseif (previous==3)&(current==1)
    transitions[7] += 1
  elseif (previous==3)&(current==2)
    transitions[8] += 1
  elseif (previous==3)&(current==-999)
    transitions[9] += 1
  else
    # do nothing.
  end
  return transitions
end




"""
`CondSum(hos::hospital; DRG = 7)`
For each DRG - need a conditional sum at each level.
times two types of patients.
"""
function CondSum(hos::hospital; DRG = 7)
  len = size(hos.levelhistory, 1)
  private = zeros(Int64, DRG,3)
  medicaid = zeros(Int64, DRG,3)
  wtp_out = zeros(Float64, DRG, 3)
  transitions = zeros(Int64, 9) # record transitions
  for el in 1:len
    if hos.levelhistory[el] == 1
      ##########  NB: Begin Private Section #######
      private[1,1] += hos.pdemandhist.demand385[el]
      private[2,1] += hos.pdemandhist.demand386[el]
      private[3,1] += hos.pdemandhist.demand387[el]
      private[4,1] += hos.pdemandhist.demand388[el]
      private[5,1] += hos.pdemandhist.demand389[el]
      private[6,1] += hos.pdemandhist.demand390[el]
      private[7,1] += hos.pdemandhist.demand391[el]
      ##########  NB: Begin Medicaid Section #######
      medicaid[1,1] += hos.mdemandhist.demand385[el]
      medicaid[2,1] += hos.mdemandhist.demand386[el]
      medicaid[3,1] += hos.mdemandhist.demand387[el]
      medicaid[4,1] += hos.mdemandhist.demand388[el]
      medicaid[5,1] += hos.mdemandhist.demand389[el]
      medicaid[6,1] += hos.mdemandhist.demand390[el]
      medicaid[7,1] += hos.mdemandhist.demand391[el]
      ########### NB: WTP Section ###########
      wtp_out[1,1] += hos.wtphist.w385[el]
      wtp_out[2,1] += hos.wtphist.w386[el]
      wtp_out[3,1] += hos.wtphist.w387[el]
      wtp_out[4,1] += hos.wtphist.w388[el]
      wtp_out[5,1] += hos.wtphist.w389[el]
      wtp_out[6,1] += hos.wtphist.w390[el]
      wtp_out[7,1] += hos.wtphist.w391[el]
    elseif hos.levelhistory[el] == 2
      ##########  NB: Begin Private Section #######
      private[1,2] += hos.pdemandhist.demand385[el]
      private[2,2] += hos.pdemandhist.demand386[el]
      private[3,2] += hos.pdemandhist.demand387[el]
      private[4,2] += hos.pdemandhist.demand388[el]
      private[5,2] += hos.pdemandhist.demand389[el]
      private[6,2] += hos.pdemandhist.demand390[el]
      private[7,2] += hos.pdemandhist.demand391[el]
      ##########  NB: Begin Medicaid Section #######
      medicaid[1,2] += hos.mdemandhist.demand385[el]
      medicaid[2,2] += hos.mdemandhist.demand386[el]
      medicaid[3,2] += hos.mdemandhist.demand387[el]
      medicaid[4,2] += hos.mdemandhist.demand388[el]
      medicaid[5,2] += hos.mdemandhist.demand389[el]
      medicaid[6,2] += hos.mdemandhist.demand390[el]
      medicaid[7,2] += hos.mdemandhist.demand391[el]
      ########### NB: WTP Section ###########
      wtp_out[1,2] += hos.wtphist.w385[el]
      wtp_out[2,2] += hos.wtphist.w386[el]
      wtp_out[3,2] += hos.wtphist.w387[el]
      wtp_out[4,2] += hos.wtphist.w388[el]
      wtp_out[5,2] += hos.wtphist.w389[el]
      wtp_out[6,2] += hos.wtphist.w390[el]
      wtp_out[7,2] += hos.wtphist.w391[el]
    elseif hos.levelhistory[el] == 3
      ##########  NB: Begin Private Section #######
      private[1,3] += hos.pdemandhist.demand385[el]
      private[2,3] += hos.pdemandhist.demand386[el]
      private[3,3] += hos.pdemandhist.demand387[el]
      private[4,3] += hos.pdemandhist.demand388[el]
      private[5,3] += hos.pdemandhist.demand389[el]
      private[6,3] += hos.pdemandhist.demand390[el]
      private[7,3] += hos.pdemandhist.demand391[el]
      ##########  NB: Begin Medicaid Section #######
      medicaid[1,3] += hos.mdemandhist.demand385[el]
      medicaid[2,3] += hos.mdemandhist.demand386[el]
      medicaid[3,3] += hos.mdemandhist.demand387[el]
      medicaid[4,3] += hos.mdemandhist.demand388[el]
      medicaid[5,3] += hos.mdemandhist.demand389[el]
      medicaid[6,3] += hos.mdemandhist.demand390[el]
      medicaid[7,3] += hos.mdemandhist.demand391[el]
      ########### NB: WTP Section ###########
      wtp_out[1,3] += hos.wtphist.w385[el]
      wtp_out[2,3] += hos.wtphist.w386[el]
      wtp_out[3,3] += hos.wtphist.w387[el]
      wtp_out[4,3] += hos.wtphist.w388[el]
      wtp_out[5,3] += hos.wtphist.w389[el]
      wtp_out[6,3] += hos.wtphist.w390[el]
      wtp_out[7,3] += hos.wtphist.w391[el]
    else # -999 - exited.
      # skip
    end
    ##########  NB: Begin Transition Record Section #######
    if el > 1
      if hos.levelhistory[el-1] != hos.levelhistory[el]
        transitions += TransitionGen(hos.levelhistory[el], hos.levelhistory[el-1])
      end
    end
  end
  # Reshape these before returning - they are now vectors with 7 entries for each DRG times 3 levels.
  # return private, medicaid, wtp_out, transitions
  return reshape(private, 1, DRG*3), reshape(medicaid, 1, DRG*3), reshape(wtp_out, 1, DRG*3), transitions'
end



"""
`DemandCheck(Tex::EntireState)`
Not so useful - just prints everyone's history of demand and the sum
"""
function DemandCheck(Tex::EntireState)
  for el in Tex.ms
    for hos in keys(el.collection)
      println(el.collection[hos].pdemandhist)
      println(CondSum(el.collection[hos]))
    end
  end
end





"""
`ResultsOut(Tex::EntireState, OtherTex::EntireState; T::Int64 = 50, beta::Float64 = 0.95,  dim2::Int64 = 81, drgamt::Array{Float64,1} = [12038.83, 66143.19, 19799.52, 4044.67, 6242.39, 1329.98, 412.04])`
Maps all of the hospital results out to a big matrix, sorted in the first column by the fid.
"""
function ResultsOut(Tex::EntireState, OtherTex::EntireState; T::Int64 = 50, beta::Float64 = 0.95,  dim2::Int64 = 81, drgamt::Array{Float64,1} = [12038.83, 66143.19, 19799.52, 4044.67, 6242.39, 1329.98, 412.04]) #dim2 - 33 paramsx2 + 7x2 records of medicaid volumes + one identifying FID
  dim1 = Tex.fipsdirectory.count
  outp = Array{Float64,2}(dim1, dim2)
  fids = [k for k in keys(Tex.fipsdirectory)]
  for el in 1:size(fids,1)
    outp[el,1] = fids[el]                                                           # Write out all of the fids as an ID in the first column.
  end
  for el in keys(Tex.fipsdirectory)                                                 # Now this is all of the hospitals by FID
    hosp = Tex.mkts[Tex.fipsdirectory[el]].collection[el]
    outprob = prod(hosp.probhistory)                                                # Prob of the outcome.
    private, medicaid, wtp_out, transitions = CondSum(hosp)
    arr = zeros(1, 40)
    arr[1] = (beta^T)*outprob*dot(wtp_out[1:7], private[1:7])                       # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 1
    arr[2] = (beta^T)*outprob*dot(wtp_out[8:14], private[8:14])                     # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 2
    arr[3] = (beta^T)*outprob*dot(wtp_out[15:21], private[15:21])                   # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 3
    arr[4] = (beta^T)*outprob*(private[1]+medicaid[1])                               # The next lines are patients summed over types.  Costs are treated as the same over Medicaid and privately insured.
    arr[5] = (beta^T)*outprob*(private[8]+medicaid[8])
    arr[6] = (beta^T)*outprob*(private[15]+medicaid[15])
    arr[7] = (beta^T)*outprob*(private[2]+medicaid[2])
    arr[8] = (beta^T)*outprob*(private[9]+medicaid[9])
    arr[9] = (beta^T)*outprob*(private[16]+medicaid[16])
    arr[10] = (beta^T)*outprob*(private[3]+medicaid[3])
    arr[11] = (beta^T)*outprob*(private[10]+medicaid[10])
    arr[12] = (beta^T)*outprob*(private[17]+medicaid[17])
    arr[13] = (beta^T)*outprob*(private[4]+medicaid[4])
    arr[14] = (beta^T)*outprob*(private[11]+medicaid[11])
    arr[15] = (beta^T)*outprob*(private[18]+medicaid[18])
    arr[16] = (beta^T)*outprob*(private[5]+medicaid[5])
    arr[17] = (beta^T)*outprob*(private[12]+medicaid[12])
    arr[18] = (beta^T)*outprob*(private[19]+medicaid[19])
    arr[19] = (beta^T)*outprob*(private[6]+medicaid[6])
    arr[20] = (beta^T)*outprob*(private[13]+medicaid[13])
    arr[21] = (beta^T)*outprob*(private[20]+medicaid[20])
    arr[22] = (beta^T)*outprob*(private[7]+medicaid[7])
    arr[23] = (beta^T)*outprob*(private[14]+medicaid[14])
    arr[24] = (beta^T)*outprob*(private[21]+medicaid[21])
    arr[25] = (beta^T)*outprob*(medicaid[1]+medicaid[8]+medicaid[15])*drgamt[1]    # Patients*revenue avg. at DRG 385
    arr[26] = (beta^T)*outprob*(medicaid[2]+medicaid[9]+medicaid[16])*drgamt[2]
    arr[27] = (beta^T)*outprob*(medicaid[3]+medicaid[10]+medicaid[17])*drgamt[3]
    arr[28] = (beta^T)*outprob*(medicaid[4]+medicaid[11]+medicaid[18])*drgamt[4]
    arr[29] = (beta^T)*outprob*(medicaid[5]+medicaid[12]+medicaid[19])*drgamt[5]
    arr[30] = (beta^T)*outprob*(medicaid[6]+medicaid[13]+medicaid[20])*drgamt[6]
    arr[31] = (beta^T)*outprob*(medicaid[7]+medicaid[14]+medicaid[21])*drgamt[7]
    arr[32:end] = (alltrans = (beta^T)*outprob*transitions)
    index = findfirst(outp[:,1], hosp.fid)                                    # find where the fid is in the list.
    outp[index, 2:41] = arr
    # NB: Here starts the second state record.
    hosp_neq = OtherTex.mkts[OtherTex.fipsdirectory[el]].collection[el]       # Find the record in the OTHER EntireState
    outprobn = prod(hosp_neq.probhistory)                                     # Prob of the outcome.
    privaten, medicaidn, wtp_outn, transitionsn = CondSum(hosp_neq)
    narr = zeros(1, 40)
    narr[1] = (beta^T)*outprobn*dot(wtp_outn[1:7], privaten[1:7])             # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 1
    narr[2] = (beta^T)*outprobn*dot(wtp_outn[8:14], privaten[8:14])           # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 2
    narr[3] = (beta^T)*outprobn*dot(wtp_outn[15:21], privaten[15:21])         # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 3
    narr[4] = (beta^T)*outprobn*(privaten[1]+medicaidn[1])                    # The next lines are patients summed over types.  Costs are treated as the same over Medicaid and privately insured.
    narr[5] = (beta^T)*outprobn*(privaten[8]+medicaidn[8])
    narr[6] = (beta^T)*outprobn*(privaten[15]+medicaidn[15])
    narr[7] = (beta^T)*outprobn*(privaten[2]+medicaidn[2])
    narr[8] = (beta^T)*outprobn*(privaten[9]+medicaidn[9])
    narr[9] = (beta^T)*outprobn*(privaten[16]+medicaidn[16])
    narr[10] = (beta^T)*outprobn*(privaten[3]+medicaidn[3])
    narr[11] = (beta^T)*outprobn*(privaten[10]+medicaidn[10])
    narr[12] = (beta^T)*outprobn*(privaten[17]+medicaidn[17])
    narr[13] = (beta^T)*outprobn*(privaten[4]+medicaidn[4])
    narr[14] = (beta^T)*outprobn*(privaten[11]+medicaidn[11])
    narr[15] = (beta^T)*outprobn*(privaten[18]+medicaidn[18])
    narr[16] = (beta^T)*outprobn*(privaten[5]+medicaidn[5])
    narr[17] = (beta^T)*outprobn*(privaten[12]+medicaidn[12])
    narr[18] = (beta^T)*outprobn*(privaten[19]+medicaidn[19])
    narr[19] = (beta^T)*outprobn*(privaten[6]+medicaidn[6])
    narr[20] = (beta^T)*outprobn*(privaten[13]+medicaidn[13])
    narr[21] = (beta^T)*outprobn*(privaten[20]+medicaidn[20])
    narr[22] = (beta^T)*outprobn*(privaten[7]+medicaidn[7])
    narr[23] = (beta^T)*outprobn*(privaten[14]+medicaidn[14])
    narr[24] = (beta^T)*outprobn*(privaten[21]+medicaidn[21])
    narr[25] = (beta^T)*outprobn*(medicaidn[1]+medicaidn[8]+medicaidn[15])*drgamt[1]
    narr[26] = (beta^T)*outprobn*(medicaidn[2]+medicaidn[9]+medicaidn[16])*drgamt[2]
    narr[27] = (beta^T)*outprobn*(medicaidn[3]+medicaidn[10]+medicaidn[17])*drgamt[3]
    narr[28] = (beta^T)*outprobn*(medicaidn[4]+medicaidn[11]+medicaidn[18])*drgamt[4]
    narr[29] = (beta^T)*outprobn*(medicaidn[5]+medicaidn[12]+medicaidn[19])*drgamt[5]
    narr[30] = (beta^T)*outprobn*(medicaidn[6]+medicaidn[13]+medicaidn[20])*drgamt[6]
    narr[31] = (beta^T)*outprobn*(medicaidn[7]+medicaidn[14]+medicaidn[21])*drgamt[7]
    narr[32:end] = (beta^T)*outprobn*transitionsn
    outp[index, 42:end] = narr
  end
  return sortrows(outp, by=x->x[1])                                                    # sort by first column (fid)
end





"""
`OuterSim(MCcount::Int; T1::Int64 = 3, dim1::Int64 = 290, dim2::Int64 = 67, fi = fips, ProjectModule.alldists)`
Runs the Monte Carlo - Equilibrium and Non-equilibrium simulations for each market MCcount times.
Note that the reduction is (+), but that includes adding the fids, so this must be divided by MCcount
to return correct results.
"""
function OuterSim(MCcount::Int; T1::Int64 = 3, dim1::Int64 = 290, dim2::Int64 = 67, fi = ProjectModule.fips, da = ProjectModule.alldists)
  outp = @sync @parallel (+) for j = 1:MCcount
    println("Current iteration ", j)
    TexasEq = MakeNew(fi, da);                                                                           # Returns an EntireState.  very quick ≈ 0.1 seconds.
    #TexasNeq = MakeNew(fi, da);                                                                         # Returns a separate EntireState.
    eq_patients = NewPatients(TexasEq)                                                                   # Separate patients - these linked to Eq Entire State.
    #neq_patients = NewPatients(TexasNeq)                                                                # Separate patients - these linked to Neq Entire State.
    ResultsOut(NewSim(T1, TexasEq, eq_patients), PSim(T1); T = T1)                                       # simulates and writes out the results.
  end
  outp[:,1] = outp[:,1]/MCcount                                                                          # Combined by (+) so reproduce the fids by dividing.
  return outp
end




"""
`ZeroFind(mat::Array{Float64,2})`
looks for hospitals which were never demanded - there are some which have rows of zeros
will want to pop those out of the results matrix.
"""
function ZeroFind(mat::Array{Float64,2})
  zers = zeros(size(mat,1)) # vector of zeros for output
  for i = 1:size(mat,1)
    zercnt = 0
    for j = 1:size(mat, 2) # don't want to count all of them - fix this.
      if mat[i,j] == 0.0
        zercnt += 1
      end
    end
    zers[i] = zercnt
  end
  return zers
end



## NB: To clean up - but these aren't debugged or used.




"""
`HospitalClean(hos::hospital)`
resets the hospital to the initial state after a run of the simulation.
some fields don't change: fid, lat, long, name, fips, bedcount,
eventually perturbed will probably change.
"""
function HospitalClean(hos::hospital)
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




"""
`Restore(Tex::EntireState)`
This function needs to set all hospital states back to zero.
Then it re-computes the set of neighbors using NeighborFix, fixing hosp.neigh and hosp.hood.
Finally it recomputes the initial choice probabilities.
This does not work yet.
"""
function Restore(Tex::EntireState)
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



###
