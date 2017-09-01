#=

works on:
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

        ##### NB: Supply-side Data Creation Functions ######




"""
`MakeIt(Tex::EntireState, fip::Vector)`

Testing this:
Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
MakeIt(Tex, ProjectModule.fips);

"""
function MakeIt(Tex::EntireState, fip::Vector{Int64})
  for f in fip
    if f != 0
      push!(Tex.ms, Market( Array{hospital,1}(), Dict{Int64, hospital}(), f, Dict{Int64, Bool}()))
    end
  end
  for m in Tex.ms # comprehensions don't work with immutables.  
    Tex.mkts[m.fipscode] = m
  end 
  return Tex
end



"""
`TXSetup(Tex::EntireState, data::Matrix, sp::Int64 ; ...)`
 Takes an entire state and adds data from the imported choices returns a record with
 fipscodes containing hospitals with mostly empty field values.  Arrays will be length sp.

 Testing:

 Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
 MakeIt(Tex, ProjectModule.fips);
 TXSetup(Tex, ProjectModule.alldists, 12);

"""
function TXSetup(Tex::EntireState,
                 data::Matrix,
                 sp::Int64;
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
                  initial(level) , #initial level - immutable.
                  Array{Int64,1}(sp), #level history
                  DemandHistory( Array{Int64,1}(sp),  Array{Int64,1}(sp), Array{Int64,1}(sp), Array{Int64,1}(sp), Array{Int64,1}(sp),  Array{Int64,1}(sp), Array{Int64,1}(sp) ), #private demand history
                  DemandHistory( Array{Int64,1}(sp),  Array{Int64,1}(sp), Array{Int64,1}(sp), Array{Int64,1}(sp), Array{Int64,1}(sp),  Array{Int64,1}(sp), Array{Int64,1}(sp) ), #medicaid demand history
                  WTP( Array{Float64,1}(sp),  Array{Float64,1}(sp), Array{Float64,1}(sp), Array{Float64,1}(sp), Array{Float64,1}(sp),  Array{Float64,1}(sp), Array{Float64,1}(sp) ),
                  Weights([0.0]), #choice probs
                  Array{Float64,1}(sp), #prob history
                  neighbors(0,0,0,0,0,0,0,0,0), #neighbors
                  Array{Int64,1}(), # hood (array of fids) - uncertain length.
                    0    , # beds added later.
                  false ) ) # perturbed or not?
      end
    end
    # push all hospital fid/ fips pairs into the directory.
    Tex.fipsdirectory[data[i, fidcol]] = fips # now for the whole state I can immediately figure out which market a hospital is in.
  end
  NeighborFix(Tex) #TODO 02/14/17 check length issues.
  InitChoice(Tex) #TODO 02/14/17 check length issues.
  ExpandDict(Tex) #TODO 02/14/17 check length issues.
  return Tex
end

"""
`FindFids(m::Market)`
Return all the fids in the market.
This is necessary in TXSetup to make sure that we only add each hospital to the market once.

Testing:
  Setup -
  #NB: Make It Needs a state argument.
  Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
  MakeIt(Tex, ProjectModule.fips);
Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists, 50);
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
#NB: Make It Needs a state argument.
Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
MakeIt(Tex, ProjectModule.fips);
Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists, 50);
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
        el.chprobability = Weights(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
      elseif el.level == 2
        levl = (1,0)
        levels = MktSize(el.neigh)
        el.chprobability = Weights(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
      elseif el.level == 3
        levl = (0,1)
        levels = MktSize(el.neigh)
        el.chprobability = Weights(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
      end
#      levels = MktSize(el.neigh)
#      el.chprobability = Weights(vec(logitest(levl, levels[1], levels[2], levels[3], [el.neigh.level105; el.neigh.level205; el.neigh.level305; el.neigh.level1515; el.neigh.level2515; el.neigh.level3515; el.neigh.level11525; el.neigh.level21525; el.neigh.level31525 ] )))
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
    # TODO 02/14/2017 - these comprehensions are surely slow.
    el.collection = Dict(i.fid => i for i in el.config)
    # for i in el.config 
    #     el.collection[i.fid] = i
    # end 

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


#=
"""
NB - this function is being deprecated.
`MakeNew(fi::Vector, dat::Matrix)`
Call this and the whole state with all markets should be created.
Should be called on "fips" or ProjectModule.fips and ProjectModule.alldists.
MakeNew(ProjectModule.fips, ProjectModule.alldists)

"""
function MakeNew(fi::Vector, dat::Matrix)
  Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
  MakeIt(Tex, fi)
  TXSetup(Tex, dat)
  return Tex
end
=#


"""
`CreateEmpty(fi::Vector, dat::Matrix,sp::Int64)`
 This creates an empty entire state record for the perturbed simulation.

Testing:
Tex = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50)

"""
function CreateEmpty(fi::Vector, dat::Matrix, sp::Int64)
  Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
  MakeIt(Tex, fi)
  TXSetup(Tex, dat, sp)
  return Tex
end

#NB: Creates hospital datastructure
#Texas = MakeNew(fips, ProjectModule.alldists );

      #### NB:  Supply-side  Printing Utilities to Display Simulation Outcomes




  ### NB: Substantive Supply-side Functions.




"""
`NewEntrantLocation(mkt::Market)`
Takes the market, takes the mean location of all hospitals, adds normal noise to it.  ≈ 6 miles perturbation from mean.
NB: return of Union{,} comes from rand(Normal(0, 0.1))
"""
function NewEntrantLocation(mkt::Market)
  meanlat::Float64 = 0.0
  meanlong::Float64 = 0.0
  m_sz::Float64 = convert(Float64, length(mkt.config))
  for el in mkt.config # over hospitals
    meanlat += el.lat
    meanlong += el.long
  end
  meanlat/=m_sz
  meanlong/=m_sz
  meanlat += rand(Normal(0, 0.1))
  meanlong += rand(Normal(0, 0.1))
  return [meanlat, meanlong]::Array{Float64,1}
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

Testing this:
#NB: Make It Needs a state argument.
Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists, 50);
Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
MakeIt(Tex, ProjectModule.fips);
NeighborRemove(Tex.mkts[48453].config[1], Tex.mkts[48453].config[2])
NeighborAppend(Tex.mkts[48453].config[1], Tex.mkts[48453].config[2])

"""
function NeighborAppend{T<:Fac}(elm::T, entrant::T)
  dist::Float64 = distance(elm.lat, elm.long, entrant.lat, entrant.long )
  if !in(entrant.fid, elm.hood)
    if (dist < 25.0)&(entrant.level != -999)
      push!(elm.hood, entrant.fid)
      if dist<5.0
        if entrant.level == 1
          elm.neigh.level105 += 1
        elseif entrant.level == 2
          elm.neigh.level205 += 1
        elseif entrant.level == 3
          elm.neigh.level305 += 1
        end
      elseif (dist>5.0)&(dist<15.0)
        if entrant.level == 1
          elm.neigh.level1515 += 1
        elseif entrant.level == 2
          elm.neigh.level2515 += 1
        elseif entrant.level == 3
          elm.neigh.level3515 += 1
        end
      elseif (dist>15.0)
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
Testing this:
#NB: Make It Needs a state argument.

Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists, 50);
Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
MakeIt(Tex, ProjectModule.fips);
NeighborRemove(Tex.mkts[48453].config[1], Tex.mkts[48453].config[2])
NeighborAppend(Tex.mkts[48453].config[1], Tex.mkts[48453].config[2])

"""
function NeighborRemove{T<:Fac}(elm::T, entrant::T)
  dist::Float64 = distance(elm.lat, elm.long, entrant.lat, entrant.long )
  if in(entrant.fid, elm.hood)
    if dist < 25.0
      deleteat!(elm.hood, findin(elm.hood, entrant.fid))
      if dist<5.0
        if entrant.level == 1
          elm.neigh.level105 = max(elm.neigh.level105 -1, 0)
        elseif entrant.level == 2
          elm.neigh.level205 = max(elm.neigh.level205 -1, 0)
        elseif entrant.level == 3
          elm.neigh.level305 = max(elm.neigh.level305 -1, 0)
        end
      elseif (dist>5.0)&(dist<15.0)
        if entrant.level == 1
          elm.neigh.level1515 = max(elm.neigh.level1515 -1,0)
        elseif entrant.level == 2
          elm.neigh.level2515 = max(elm.neigh.level2515 -1,0)
        elseif entrant.level == 3
          elm.neigh.level3515 = max(elm.neigh.level3515 -1,0)
        end
      elseif (dist>15.0)
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
      # avoid allocations.
      hosp.neigh.level105 = 0
      hosp.neigh.level205 = 0
      hosp.neigh.level305 = 0
      hosp.neigh.level1515 = 0
      hosp.neigh.level2515 = 0
      hosp.neigh.level3515 = 0
      hosp.neigh.level11525 = 0
      hosp.neigh.level21525 = 0
      hosp.neigh.level31525 = 0
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
  for mk1 in state.ms #over market fips
    for h1 in mk1.config #hospitals in the market
      for mk2 in state.ms # now markets in the state again.
        for h2 in mk2.config # hospitals in the market.
          if (h1.fid != h2.fid)&&(distance(h1.lat, h1.long, h2.lat, h2.long) < 25.0) # if they are closer than 25 AND not the same facility.
            NeighborAppend(h1, h2)
            NeighborAppend(h2, h1)
          end
        end
      end
    end
  end
end






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
    delete!(mkt.collection, el)
  end
end




"""
`HospUpdate(hosp::hospital, choice::Int; update = false)`
Takes a hospital record and updates the probabilities of the choices.

Testing:
#Tex = TXSetup(MakeIt(ProjectModule.fips), ProjectModule.alldists, 50);
Tex = EntireState(Array{hospital,1}(), Dict{Int64,Market}(), Dict{Int64,hospital}())
MakeIt(Tex, ProjectModule.fips);
TXSetup(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);
NewSim(10, Tex, patients);

HospUpdate(Texas.mkts[48453].config[1], 1; update = true)
"""
function HospUpdate{T<:ProjectModule.Fac}(hosp::T, choice::Int64; update = false)
 levl = (-1, -1)
 if (hosp.level!=choice)||update # want to be able to force this to rerun when the data is cleaned again.
   if choice != -999
     if choice == 1
       levl = (0,0)
     elseif choice == 2
       levl = (1,0)
     elseif choice == 3
       levl = (0,1)
     end
     # FIXME - error here.  What is output of logitest?  What is type of input to Weights?  
     levels = MktSize(hosp.neigh)
     prs = logitest(levl, levels[1], levels[2], levels[3], [hosp.neigh.level105; hosp.neigh.level205; hosp.neigh.level305; hosp.neigh.level1515; hosp.neigh.level2515; hosp.neigh.level3515; hosp.neigh.level11525; hosp.neigh.level21525; hosp.neigh.level31525 ] )
     return Weights(vec(prs))
   else # choice = -999
     return Weights([1.0]) #TODO: one option, no choices ??  Might need four options [1.0 1.0 1.0 1.0]
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
     return Weights(vec(prs))
   else # choice = -999
     return Weights([1.0]) #TODO: one option, no choices ??  Might need four options [1.0 1.0 1.0 1.0]
   end
 else
   return hosp.chprobability
  end
end

      ##### NB: Demand-side Data Structure Creation.


"""
`function CreateZips(alld::Array,Tex::EntireState)::patientcollection`
This creates zip code records but doesn't put patients into them.
Those tasks are done by FillMPatients and FillPPatients.
This contains zips which are *not* seen in the data.  That's potentially ok, but also potentially confusing.
But check to make sure *some* of the zips have patients.
Now the correct distances are being imported as ProjectModule.alldists.
But in this version do check to make sure that the distance between any fac and zip is
less than 25 miles.
Note that this indexes into the Array according to the column zipcol.

Testing:
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
zipcheck = CreateZips(ProjectModule.alldists, Tex);

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
  # TODO - hard-code the model coefficients here.  
  ppatients::patientcollection = patientcollection( Dict{Int64, zipcode}() )
  # unfound = Array{Int64,1}()
  # found = Array{Int64,1}()
  # names = Array{AbstractString, 1}()
  fids = unique(dat[:,fidcol])
  matched::Int64 = 0
  # TODO - this comprehension is slow.  
  ppatients.zips = Dict(k=> zipcode(k, 0, Dict{Int64,ProjectModule.Fac}(), Dict{Int64,Float64}(),                                                                              # zipcode, public health region, facilities, hospital FE's.
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
`FillPatients(pats::patientcollection, pcounts::Matrix)`
Adds the privately insured and medicaid patients to the zip records.
Combines in one function the previous FillPPatients and FillMPatients.
Uses different data too - pcounts rather than Medicaid/Private Individual Choices.  

Testing: 

Texas2 = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients2 = CreateZips(ProjectModule.alldists, Texas2);
FillPatients(patients2, ProjectModule.pcount)


"""
function FillPatients(pats::patientcollection, pcounts::Matrix)
  # use 2005 medicaid and private numbers.  These are columns 19 (private), 20 (medicaid)
  mloc::Int64 = 20
  ploc::Int64 = 19
  ziploc::Int64 = 1
  drgloc::Int64 = 2
  notavail::Array{Float64,1} = setdiff(ProjectModule.pcount[:,ziploc],keys(pats.zips))
  for r in 1:size(pcounts,1)
    if !in(pcounts[r,ziploc], notavail)
      if pcounts[r,drgloc] == 385
        pats.zips[pcounts[r,ziploc]].mpatients.count385 = pcounts[r,mloc] 
        pats.zips[pcounts[r,ziploc]].ppatients.count385 = pcounts[r,ploc]
      elseif pcounts[r,drgloc] == 386 
        pats.zips[pcounts[r,ziploc]].mpatients.count386 = pcounts[r,mloc] 
        pats.zips[pcounts[r,ziploc]].ppatients.count386 = pcounts[r,ploc]
      elseif pcounts[r,drgloc] == 387 
        pats.zips[pcounts[r,ziploc]].mpatients.count387 = pcounts[r,mloc] 
        pats.zips[pcounts[r,ziploc]].ppatients.count387 = pcounts[r,ploc]
      elseif pcounts[r,drgloc] == 388 
        pats.zips[pcounts[r,ziploc]].mpatients.count388 = pcounts[r,mloc] 
        pats.zips[pcounts[r,ziploc]].ppatients.count388 = pcounts[r,ploc]
      elseif pcounts[r,drgloc] == 389
        pats.zips[pcounts[r,ziploc]].mpatients.count389 = pcounts[r,mloc] 
        pats.zips[pcounts[r,ziploc]].ppatients.count389 = pcounts[r,ploc]
      elseif pcounts[r,drgloc] == 390
        pats.zips[pcounts[r,ziploc]].mpatients.count390 = pcounts[r,mloc] 
        pats.zips[pcounts[r,ziploc]].ppatients.count390 = pcounts[r,ploc]
      elseif pcounts[r,drgloc] == 391
        pats.zips[pcounts[r,ziploc]].mpatients.count391 = pcounts[r,mloc] 
        pats.zips[pcounts[r,ziploc]].ppatients.count391 = pcounts[r,ploc]
      else 
        println("fail - ", pcounts[r, drgloc], " ", pcounts[r, ziploc])
      end 
    end 
  end 
end 




"""
`ComputeDetUtil(zipc::zip, fid::Int64, p_or_m::Bool)`
Computes the deterministic component of utility for each hospital in the zip "zipc".
Maps exited facilites to have deterministic utility -999
Works on private and medicaid patients by setting p_or_m to true or false, respectively.
Has been written to accomodate hospital FE's when available.
Allocations come here from accessing the fields.  18-20 per call.
0.000014 seconds (18 allocations: 384 bytes) -> this is one call in one of the lines below, not the whole thing.
"""
function ComputeDetUtil(zipc::zipcode, fid::Int64, p_or_m::Bool)
  dist = distance(zipc.facilities[fid].lat, zipc.facilities[fid].long, zipc.lat, zipc.long)
  if p_or_m #if TRUE private
    if zipc.facilities[fid].level == 1
      zipc.pdetutils[fid] = zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)#+zipc.pcoeffs.closest*(0) + zipc.fes[fid]
    elseif zipc.facilities[fid].level == 2
      zipc.pdetutils[fid] = zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.inter+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)#+zipc.pcoeffs.closest*(0) + zipc.fes[fid]
    elseif zipc.facilities[fid].level == 3
      zipc.pdetutils[fid] = zipc.pcoeffs.distance*dist+zipc.pcoeffs.distsq*(dist^2)+zipc.pcoeffs.inten+zipc.pcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)#+zipc.pcoeffs.closest*(0) + zipc.fes[fid]
    else  # =-999
      zipc.pdetutils[fid] = -999.0 # can't choose a facility which has exited - set det utility very low.
    end
  else
    if zipc.facilities[fid].level == 1
      zipc.mdetutils[fid] = zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)#+zipc.pcoeffs.closest*(0) + zipc.fes[fid]
    elseif zipc.facilities[fid].level == 2
      zipc.mdetutils[fid] = zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.inter+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)#+zipc.pcoeffs.closest*(0) + zipc.fes[fid]
    elseif zipc.facilities[fid].level == 3
      zipc.mdetutils[fid] = zipc.mcoeffs.distance*dist+zipc.mcoeffs.distsq*(dist^2)+zipc.mcoeffs.inten+zipc.mcoeffs.distbed*(dist*zipc.facilities[fid].bedcount/100)#+zipc.pcoeffs.closest*(0) + zipc.fes[fid]
    else  # =-999
      zipc.mdetutils[fid] = -999.0
    end
  end
end





"""
`UpdateDeterministic(collect::patientcollection)`
Computes the deterministic component of the utility - updates every firm every time it is called.
Is called during the Eq and Non-eq simulations.
0.043695 seconds (483.33 k allocations: 7.605 MB)

## Testing ##

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
UpdateDeterministic(patients)

"""
function UpdateDeterministic(collt::patientcollection)
  for el in keys(collt.zips) #iterates over zips
    for fid in keys(collt.zips[el].facilities) # iterates over dict of facilities within zip.
      ComputeDetUtil(collt.zips[el], fid, false)
      ComputeDetUtil(collt.zips[el], fid, true)
    end
  end
end




"""
`NewPatients(Tex::EntireState; dists = ProjectModule.alldists, phrloc = 103, pins = pinsured, pmed = pmedicaid)`
this creates the whole collection of patients.
It must take an existing EntireState record to link the hospitals.

#TODO - there are patients who don't get added.  Their zips are not in the list but are in the inpatient discharge.

Testing:
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);

Timing:
0.248355 seconds (2.75 M allocations: 54.025 MB)
"""
function NewPatients(Tex::EntireState;
                     dists::Array{Any,2} = ProjectModule.alldists,
                     pins = ProjectModule.pcount)
  patients = CreateZips(dists, Tex) #NB: This needs to take the whole state so that the hosps in zips point to the same underlying record.
  FillPatients(patients, pins)
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
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
NewPatientsTest(Texas)

# Results - a bunch of missing zips.  Try to load and find them in that third source of zip codes.
# But how many patients are there in these zips?  That is not clear yet.  Figure that out first.
And what happens to them in the previous function when the set of patients is created?  Not sure.

Do I have locations for all of the zips?  That is the question.
- Check CreateZips.  Around line 700.
"""
function NewPatientsTest(Tex::EntireState;
                     dists::Array{Any,2} = ProjectModule.alldists,
                     phrloc::Int64 = 103,
                     ziploc::Int64 = 101,
                     pins = ProjectModule.pcount)
  patients = CreateZips(dists, Tex) # creates empty zip codes.
  patients = FillPatients(patients, pins)
  UpdateDeterministic(patients)
  AddOO(patients)
  pzips::Dict{Int64, Int64} = Dict(k => 0 for k in unique(pins[:,ziploc])) # this is guaranteed to get all of the zips.
  for n in 1:size(pins, 1)
    if haskey(pzips, pins[n,ziploc])
      pzips[pins[n,ziploc]] += 1
    end
  end
  mzips::Dict{Int64, Int64} = Dict(k => 0 for k in unique(pmed[:,ziploc]))
  for n in 1:size(pmed,1)
    if haskey(mzips, pmed[n,ziploc])
      mzips[pmed[n,ziploc]] += 1
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
  return unique(nothere), pzips, mzips
end


"""
`FillState(Tex::EntireState, data::Matrix, ns::Int64;fidcol::Int64 = 4,latcol::Int64 = 21,longcol::Int64 = 22,namecol::Int64 = 5,fipscol::Int64 = 7,intensivecol::Int64 = 11,intermediatecol::Int64 = 20))`
fills the entire state record with elements of the chospital type
for the counterfactual only.
Note that this needs to be called AFTER the function `CMakeIt(Tex::EntireState, fip::Vector)` is called
on an empty state record.
This version also adds all of the `chospitals` directly to the `Market.collection` dictionary.
Tex = EntireState(Array{chospital,1}(), Dict{Int64,Market}(), Dict{Int64,chospital}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
"""
function FillState(Tex::EntireState, data::Matrix, ns::Int64;
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
                  Array{Int64,1}(ns), #volume
                  Array{Int64, 1}(ns), #mortality
                  Array{Float64,1}(ns), #ppayoff
                  Array{Float64,1}(ns), #mpayoff
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
end



"""
`CMakeIt(Tex::EntireState, fip::Vector)`
Perhaps poor practice to use Eval in this way, but generates markets named m*fipscode* for any fipscode in the vector fip.
Tex = EntireState(Array{chospital,1}(), Dict{Int64,Market}(), Dict{Int64,chospital}())
CMakeIt(Tex, ProjectModule.fips);


"""
function CMakeIt(Tex::EntireState, fip::Vector)
  for el in fip
    if el != 0
      push!(Tex.ms, Market( Array{chospital,1}(), Dict{Int64, chospital}(), el, Dict{Int64, Bool}()))
    end
  end
  Tex.mkts = Dict(m.fipscode => m for m in Tex.ms)
end

"""
`CounterClean(Tex::EntireState)`
Takes the initial counterfactual state and writes all values of all hosp quantities to 0.
"""
function CounterClean(Tex::EntireState)
    for fips in keys(Tex.mkts)
        for fid in keys(Tex.mkts[fips].collection)
            for i = 1:length(Tex.mkts[fips].collection[fid].totalv)
                Tex.mkts[fips].collection[fid].totalv[i] = 0
                Tex.mkts[fips].collection[fid].mortality[i] = 0
                Tex.mkts[fips].collection[fid].ppayoff[i] = 0.0
                Tex.mkts[fips].collection[fid].mpayoff[i] = 0.0
            end 
        end 
    end 
end 

"""
`CounterCleanResults(ch::counterhistory)`
Arrays are allocated so set values to zeros.  
"""
function CounterCleanResults(ch::counterhistory)
    for k1 in keys(ch.hist)
        for k2 in keys(ch.hist[k1].values)
            for k3 in keys(ch.hist[k1].values[k2].hosprecord)
                for i = 1:length(ch.hist[k1].values[k2].hosprecord[k3].totbr)
                    ch.hist[k1].values[k2].hosprecord[k3].totbr[i] = 0
                    ch.hist[k1].values[k2].hosprecord[k3].totlbw[i] = 0
                    ch.hist[k1].values[k2].hosprecord[k3].totvlbw[i] = 0
                    ch.hist[k1].values[k2].hosprecord[k3].deaths[i] = 0
                    ch.hist[k1].values[k2].hosprecord[k3].profit[i] = 0.0
                end         
            end 
        end
    end 
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






  ### NB: Substantive Demand-side Functions.




"""
`WhichZips(pats::patientcollection, fid::Int64)`
Takes a patientcollection and tells me which zips have the hospital fid
"""
function WhichZips(pats::patientcollection, fid::Int64)
  for zi in keys(pats.zips)
    # TODO 02/14/2017 - replace this with haskey()
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
`CalcWTP(zipc::zip, arr::Array{Float64,2})`
Takes the deterministic component of utility for the privately insured patients and returns a WTP measure.
Output is sent to WTPMap.

Testing:
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
@benchmark CalcWTP(patients.zips[78702]) # add dollar sign before patients.

Nice speed up - from 1.5 μs to 700 ns. 

TO: 
BenchmarkTools.Trial:
  memory estimate:  16 bytes
  allocs estimate:  1
  --------------
  minimum time:     708.850 ns (0.00% GC)
  median time:      711.886 ns (0.00% GC)
  mean time:        728.721 ns (0.61% GC)
  maximum time:     45.392 μs (97.25% GC)
  --------------
  samples:          10000
  evals/sample:     140

FROM: 

@benchmark CalcWTP(patients.zips[78759])
BenchmarkTools.Trial:
  memory estimate:  1.97 KiB
  allocs estimate:  7
  --------------
  minimum time:     1.552 μs (0.00% GC)
  median time:      1.627 μs (0.00% GC)
  mean time:        1.814 μs (7.34% GC)
  maximum time:     228.033 μs (94.43% GC)
  --------------
  samples:          10000
  evals/sample:     10 

How is this actually used again?  In NewSim...

WriteWTP(WTPMap(pats, Tex), Tex, i) # and WTPMap calls CalcWTP

newarr = zeros(2,12)
CalcWTP(patients.zips[78759], newarr)
"""
function CalcWTP(zipc::zipcode, arr::Array{Float64,2})  
  interim::Float64 = 0.0
  for (i,el) in enumerate(keys(zipc.pdetutils))
    arr[1,i] = el 
    arr[2,i] = exp(zipc.pdetutils[el])
    interim += exp(zipc.pdetutils[el]) 
  end
  for k in 1:size(arr,2)
    arr[2,k] = arr[2,k]/interim
  end 
end


"""
`WTPDict(Tex::EntireState)`
Makes a dict of FIDS from those in the state.

## Testing ## 
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
wtd = WTPDict(Texas)

"""
function WTPDict(Tex::EntireState)
  outp::Dict{Int64,Float64} = Dict{Int64,Float64}() 
  for k1 in keys(Tex.fipsdirectory) # need a list of all fids to initialize the dictionary.
    outp[k1] = 0.0
  end 
  return outp 
end 

"""
`CleanWTPDict(d::Dict{Int64,Float64})`
Takes the WTP Dict and sets all values back to zero.  

## Testing ## 
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
wtd = WTPDict(Texas)
CleanWTPDict(wtd)


"""
function CleanWTPDict(d::Dict{Int64,Float64})
  for k1 in keys(d)
    d[k1] = 0.0
  end 
end 


"""
`WTPMap(pats::patientcollection, Tex::EntireState)`
Takes a patient collection and an entire state and returns a dict{fid, WTP}
computed by calling CalcWTP.  Right now it ignores Inf and NaN.
Input is from CalcWTP.  Output is sent to WriteWTP

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
WTPMap(patients, Texas)
FROM:
memory estimate:  2.50 MiB
allocs estimate:  10737 
mean time:        6.060 ms (3.14% GC) 

TO: 
memory estimate:  30.03 KiB
allocs estimate:  1922
mean time:        3.134 ms (0.21% GC)


wtd = WTPDict(Texas)
newarr = zeros(2,12)
WTPMap(patients, Texas, wtd, newarr)

"""
function WTPMap(pats::patientcollection, Tex::EntireState, wtpd::Dict{Int64,Float64}, arr::Array{Float64,2})   
  for zipc in keys(pats.zips)
    ArrayZero(arr)                                    # call this before CalcWTP is called on each zip code   
    CalcWTP(pats.zips[zipc], arr)                     # calculate WTP and store result in arr.
    for ix1 in 1:size(arr,2)                          # all of the fids are in the first row 
      if arr[1,ix1] != 0                              # don't try to map the OO.
        if (arr[2,ix1]!=1)&(!isnan(arr[2,ix1]))       # there should be no way for this to be 1 anyway.
          wtpd[arr[1,ix1]] += log(1/(1-arr[2,ix1]))   # add WTP to dict 
        elseif (arr[2,ix1] == 1)||(isnan(arr[2,ix1])) # if an error occurs, print it.  
            println(zipc, "  ", arr[2,ix1])
        end
      end 
    end
  end                                                 # return nothing. 
end




"""
`WriteWTP(reslt::Dict{Int64, Float64}, Tex::EntireState)`
Takes a dict of {fid, WTP} and writes it out by DRG.
Works on the output of WTPMap

## Testing ## 

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
d1 = WTPMap(patients, Texas)
WriteWTP(d1, Texas, 1)

"""
function WriteWTP(reslt::Dict{Int64, Float64}, Tex::EntireState, index::Int64)
  for els in keys(reslt)
    Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w385[index]=reslt[els]
    Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w386[index]=reslt[els]
    Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w387[index]=reslt[els]
    Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w388[index]=reslt[els]
    Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w389[index]=reslt[els]
    Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w390[index]=reslt[els]
    Tex.mkts[Tex.fipsdirectory[els]].collection[els].wtphist.w391[index]=reslt[els]
  end
end


"""
`function UMap(utils::Array{Float64,1}, fids::Array{Int64,1}, temparr::Array{Float64,1})`

This function takes:
- Array of Utils
- Array of Fids
- Temporary Array

Computes utility + random component, maps out corresponding FID.

testing:
ut = [0.1, 0.2, 0.3, 0.4, 0.5];
fi = [111, 222, 333, 444, 19];
UMap(ut, fi)

10% speedup  
"""
function UMap(utils::Array{Float64,1},fids::Array{Int64,1})::Int64
  const dist_μ::Int64 = 0
  const dist_σ::Int64 = 1
  const dist_ξ::Int64 = 0
  d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
  maxm::Int64 = 1                   # this is an index to a fid.  Will record max index.
  maxu::Float64 = 0.0               # holds the value of max utility 
  tem::Float64 = 0.0                # holds interim utility value
  for i = 1:size(utils,1)
    tem = utils[i]+rand(d)
    if tem>maxu  
        maxm = i                    # replace index if greater
        maxu = tem                  # replace max util 
    end 
  end 
  return fids[maxm]::Int64          # return fid of max util value. 
end 

"""
`DV(d::Dict{Int64, Float64})::Tuple{Array{Int64,1},Array{Float64,1}}`
This is a more efficient version of `DicttoVec`.  Takes a dictionary of {Int64,Float64}
and returns two vectors.


#Testing on Choice Data:
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
DV(patients.zips[78702].pdetutils)

@code_warntype DV(patients.zips[78702].pdetutils)

This is not type stable, which is weird...
"""
function DV(d::Dict{Int64, Float64})::Tuple{Array{Int64,1},Array{Float64,1}}
  out1::Array{Int64,1} = zeros(Int64, d.count)     #for the keys/FIDs
  out2::Array{Float64,1} = zeros(Float64, d.count) #for the utils
  for (i,k) in enumerate(keys(d))
    out1[i]::Int64 = k
    out2[i]::Float64 = d[k]
  end
  return out1, out2
end


"""
`function NewHospDict(Tex::EntireState)`
Associates with each FID a zero patientcount.
The purpose is to make a dictionary to hold the output of ChoiceVector.
"""
function NewHospDict(Tex::EntireState)
  dt::Dict{Int64, patientcount} = Dict()
  for el in keys(Tex.fipsdirectory)
    dt[el] = patientcount(0,0,0,0,0,0,0)
  end
  dt[0] = patientcount(0,0,0,0,0,0,0)
  return dt
end


"""
`function PatientsClean(dt::Dict{Int64, patientcount})`
resets the values of all patientcounts to 0.  This allows re-use of
the dictionary.  Maybe the allocation could be even lower by manually
setting all fields back to 0.  But whatever.
"""
function PatientsClean(dt::Dict{Int64, patientcount})
  for el in keys(dt)
    dt[el] = patientcount(0,0,0,0,0,0,0)
  end
end



"""
`function ChoiceVector(pd::Dict{Int64, Float64},dt::Dict{Int64, patientcount},ch::Array{Int64,1},x::patientcount)`
This should operate in-place on the dictionary dt.
The dictionary has an entry for every hospital.
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
dic1 = NewHospDict(Texas);
fids1, utils1 = DV(patients.zips[78759].pdetutils)
inpt = ones(Int64, 1550); # largest group is 1511
ChoiceVector(patients.zips[78759].pdetutils, dic1, inpt, patients.zips[78759].ppatients)

FIXME - there is a type instability here from the fact that DV?
REMEMBER TO TURN ON THREADING.



"""
function ChoiceVector(pd::Dict{Int64, Float64},
                      dt::Dict{Int64, patientcount},
                      ch::Array{Int64,1},
                      x::patientcount)
  fids::Array{Int64,1}, utils::Array{Float64,1} = DV(pd) # very quick ≈ 300 ns.
  # type instability in patientcount, it seems. 
  # also in DV.   
  for (loc, nm) in enumerate(x)
    UseThreads(ch, fids, utils, nm)                      # ≈ 179 μs, for nm = 300, ≈ 12.504 μs for nm = 20
    if loc == 1
      for i = 1:nm # looping over all elements of ch which have choices recorded, rather than over the whole vector.
        dt[ch[i]].count385 += 1 # memory estimate:  48 bytes / allocs estimate:  3
      end
    elseif loc==2
      for i = 1:nm
        dt[ch[i]].count386 += 1
      end
    elseif loc==3
      for i = 1:nm
        dt[ch[i]].count387 += 1
      end
    elseif loc==4
      for i = 1:nm
        dt[ch[i]].count388 += 1
      end
    elseif loc==5
      for i = 1:nm
        dt[ch[i]].count389 += 1
      end
    elseif loc==6
      for i = 1:nm
        dt[ch[i]].count390 += 1
      end
    elseif loc==7
      for i = 1:nm
        dt[ch[i]].count391 += 1
      end
    end
    ResVec(ch) #reset the vector. - 471.337 ns (0.00% GC)
  end
end

"""
`function UseThreads(inpt::Array{Int64,1},fids::Array{Int64,1}, utils::Array{Float64,1}, temparry::Array{Float64, 1})`
The sole purpose of this is to avoid #15276, which generates an ambiguity in the type of the arrays in `ChoiceVector`.
There is no ambiguity in this one.
https://github.com/JuliaLang/julia/issues/15276
Workaround, maybe?
See this:
https://github.com/yuyichao/explore/blob/8d52fb6caa745a658f2c9bbffd3b0f0fe4a2cc48/julia/issue-17395/scale.jl#L21

# Test. # 

To get UseThreads:
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
inpt = ones(Int64, 1550); # largest group is 1511
fids1, utils1 = DV(patients.zips[78759].pdetutils)
tar = zeros(utils1)
# nm is patients.zips[78759].ppatients.count391
UseThreads(inpt, fids1, utils1,  2)


@benchmark UseThreads(inpt, fids1, utils1, 297) # put dollar signs in front of those.  
BenchmarkTools.Trial:
  memory estimate:  48 bytes
  allocs estimate:  1
  --------------
  minimum time:     177.711 μs (0.00% GC)
  median time:      179.942 μs (0.00% GC)
  mean time:        186.356 μs (0.00% GC)
  maximum time:     514.353 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
"""
function UseThreads(inpt::Array{Int64,1},fids::Array{Int64,1},utils::Array{Float64,1}, x::Int64)
  Threads.@threads for i = 1:x
    inpt[i] = UMap(utils, fids)
  end
end

"""
`function ResVec(v::Array{Int64,1})`
This function resets the vector to its original state.  This enables it to be
re-used so it does not require so many allocations.

@benchmark ResVec(inpt) # sub dollar sign.  
BenchmarkTools.Trial:
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     470.694 ns (0.00% GC)
  median time:      471.337 ns (0.00% GC)
  mean time:        487.366 ns (0.00% GC)
  maximum time:     1.271 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     196
"""
function ResVec(v::Array{Int64, 1})
  for i = 1:length(v)
    v[i] = 1
  end
end


"""
`function GenPChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})`
Returns a dictionary of private demands for the whole state.
Takes as input a patientcollection, a dict{Int64, patientcount} and a re-usable array v.

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
dic1 = NewHospDict(Texas);
inpt = ones(Int64, 1550); # largest group is 1511
GenPChoices(patients, dic1, inpt)
"""
function GenPChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})
  for k in keys(p.zips)
    ChoiceVector(p.zips[k].pdetutils, d, v, p.zips[k].ppatients)
  end
end

"""
`function GenMChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})`
Returns a dictionary of private demands for the whole state.
Takes as input a patientcollection, a dict{Int64, patientcount} and a re-usable array v.
"""
function GenMChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})
  for k in keys(p.zips)
    ChoiceVector(p.zips[k].mdetutils, d, v, p.zips[k].mpatients)
  end
end



"""
`function DicttoVec`
Takes the dictionary of utilities and returns a vector - top row is FID's and bottom row is
utilities.  Lower allocation than doing this as a comprehension.
"""
function DicttoVec(d::Dict{Int64, Float64})
  outp::Array{Float64,2} = zeros(2, d.count)
  col::Int64 = 1
  for el in keys(d)
    outp[1, col] = el
    outp[2, col] = d[el]
    col += 1
  end
  return outp
end


# FIXME - maybe delete this?  Not clear it's faster.
function PatientCountOut(Tex::EntireState)
  outp::Dict{Int64, patientcount} = Dict{Int64, patientcount}()
  for el in keys(Tex.fipsdirectory)
    outp[el] = patientcount(0,0,0,0,0,0,0)
  end
  return outp
end



"""
`PHistoryAdd(hos::hospital, cnt::patientcount)`
Maps patientcount to the private demand history
"""
function PHistoryAdd(hos::hospital, cnt::patientcount, index::Int64)
  hos.pdemandhist.demand385[index]=cnt.count385
  hos.pdemandhist.demand386[index]=cnt.count386
  hos.pdemandhist.demand387[index]=cnt.count387
  hos.pdemandhist.demand388[index]=cnt.count388
  hos.pdemandhist.demand389[index]=cnt.count389
  hos.pdemandhist.demand390[index]=cnt.count390
  hos.pdemandhist.demand391[index]=cnt.count391
end




"""
`MHistoryAdd(hos::hospital, cnt::patientcount)`
Maps patientcount to the Medicaid demand history.
"""
function MHistoryAdd(hos::hospital, cnt::patientcount, index::Int64)
  hos.mdemandhist.demand385[index]=cnt.count385
  hos.mdemandhist.demand386[index]=cnt.count386
  hos.mdemandhist.demand387[index]=cnt.count387
  hos.mdemandhist.demand388[index]=cnt.count388
  hos.mdemandhist.demand389[index]=cnt.count389
  hos.mdemandhist.demand390[index]=cnt.count390
  hos.mdemandhist.demand391[index]=cnt.count391
end




"""
`PDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)`
Maps Private patient demand out to the state record.
Note - now cleans out the dictionary from GenPChoices.

### Testing ### 

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
dic1 = NewHospDict(Texas);
inpt = ones(Int64, 1550);
GenPChoices(patients, dic1, inpt)
PDemandMap(dic1, Texas, 1)

"""
function PDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState, index::Int64)
  for el in keys(patd)
    if el != 0 # don't map out the outside option.
      PHistoryAdd(Tex.mkts[Tex.fipsdirectory[el]].collection[el], patd[el], index)
    end
  end
  PatientsClean(patd)
end




"""
`MDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState)`
Maps Medicaid Patient demand out to the state record.
Note - now cleans out the dictionary from GenMChoices.

## Testing ##

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
dic1 = NewHospDict(Texas);
inpt = ones(Int64, 1550);
GenMChoices(patients, dic1, inpt)
MDemandMap(dic1, Texas, 1)

"""
function MDemandMap(patd::Dict{Int64, patientcount}, Tex::EntireState, index::Int64)
  for el in keys(patd)
      if el != 0 # don't map out the outside option.
        MHistoryAdd(Tex.mkts[Tex.fipsdirectory[el]].collection[el], patd[el], index)
      end
  end
  PatientsClean(patd)
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
                     weightpr = weightprobs[:,2], # TODO - these are vectors of fixed size.  Hard code them and drop the import 
                     admitprobs = nicuprobs[:,2], # TODO - also a fixed size vector.  Hard code.  
                     w1 = Weights([1-admitprobs[1],admitprobs[1]]), # next few lines unused.
                     w2 = Weights([1-admitprobs[2],admitprobs[2]]),
                     w3 = Weights([1-admitprobs[3],admitprobs[3]]),
                     w4 = Weights([1-admitprobs[4],admitprobs[4]]),
                     w5 = Weights([1-admitprobs[5],admitprobs[5]]),
                     w6 = Weights([1-admitprobs[6],admitprobs[6]]),
                     w7 = Weights([1-admitprobs[7],admitprobs[7]]), # here and below all used.
                     w8 = Weights([1-admitprobs[8],admitprobs[8]]),
                     w9 = Weights([1-admitprobs[9],admitprobs[9]]),
                     w10 = Weights([1-admitprobs[10],admitprobs[10]]),
                     w11 = Weights([1-admitprobs[11],admitprobs[11]]),
                     w12 = Weights([1-admitprobs[12],admitprobs[12]]),
                     w13 = Weights([1-admitprobs[13],admitprobs[13]]))
  outp::Dict{Int64, LBW} = Dict{Int64, LBW}()  # = Dict( k=>LBW(0,0,0,0,0,0) for k in keys(Tex.fipsdirectory)) # empty dictionary of fids/LBW record types
  for k1 in keys(ppat)
    outp[k1] = LBW(0,0,0,0,0,0)
  end 
  for el in keys(ppat)
    if el!=0
        totl = sum(ppat[el] + mpat[el])
        for i = 1:totl
          bwt = sample(bins, Weights(weightpr)) # sample from the distribution of birthweights.
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
  end
  return outp, ppat, mpat
end




"""
`AllMortality(d::Dict{Int64, LBW})`
The return of `PatientDraw` is a dictionary of {fid, LBW}.  Take that volume and convert it to a mortality
rate.  Then apply the mortality rate to the LBW record.  The elements keys(d) will be fids.
"""
function AllMortality(d::Dict{Int64, LBW}, Tex::EntireState)
  #TODO 02/14/2017 - fix this comprehension, which is surely slow.
  outp::Dict{Int64, Int64} = Dict{Int64,Int64}() # = Dict(j => 0 for j in keys(d))
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

## Testing: ## 
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 20);
patients = NewPatients(Texas);
NewSim(20, Texas, patients); # 3.386766 seconds (9.28 M allocations: 183.511 MiB, 2.48% gc time)

Sources of Allocations: 
GenPChoices: memory estimate:  9.58 MiB   allocs estimate:  511404   mean time:   143.855 ms 
UpdateDeterministic: memory estimate:  7.15 MiB allocs estimate:  453242 mean time:   22.707 ms
WriteWTP: memory estimate:  62.56 KiB allocs estimate:  4004 mean time:   366.285 μs
WTPMap: memory estimate:  30.03 KiB allocs estimate:  1922 mean time:     3.061 ms (0.19% GC) 
EntryProcess: memory estimate:  422 bytes allocs estimate:  3 mean time:   84.371 μs (0.16% GC)
CleanWTPDict: memory estimate:  0 bytes allocs estimate:  0 mean time:        8.756 μs
"""
function NewSim(T::Int, Tex::EntireState, pats::patientcollection)
  const entrants::Array{Float64,1} = [0, 1, 2, 3] 
  const entryprobs::Array{Float64,1} = [0.9895, 0.008, 0.0005, 0.002]
  d1 = NewHospDict(Tex)                                                                        # creates a dict for GenP below.
  d2 = NewHospDict(Tex)                                                                        # creates a dict for GenM below
  wtpd1 = WTPDict(Tex)                                                                         # creates a dict for WTPMap
  wtparr1 = zeros(2,12)                                                                        # temporary array for CalcWTP
  arry1 = zeros(Int64, 1550)                                                                   # allocates an array for use in GenP.  Can be re-used.
  arry2 = zeros(Int64, 1550)                                                                   # allocates an array for use in GenM.  Can be re-used.
  for i = 1:T
    # TODO - check that WTP gets updated in this function, just in case.
    WTPMap(pats, Tex, wtpd1, wtparr1)
    WriteWTP(wtpd1, Tex, i) 
    GenPChoices(pats, d1, arry1)                                                               # this now modifies the dictionary in-place
    PDemandMap(d1, Tex, i)                                                                     # and this now cleans the dictionary up at the end, setting all demands to 0.
    GenMChoices(pats, d2, arry2)                                                               # this now modifies the dictionary in-place
    MDemandMap(d2, Tex, i)                                                                     # and this now cleans the dictionary up at the end, setting all demands to 0.
    for el in Tex.ms
      EntryProcess(el, i, T)                                                                   # call the Entry Process on every market - el.   
      if el.fipscode == 48453
        println("     ", i, "   ")
      end 
      for elm in el.config
        action = StatsBase.sample( ChoicesAvailable(elm), elm.chprobability )                  # Take the action
        # TODO - here add update of WTP given an action not equal to 10.
        # FIXME 
        elm.probhistory[i] = elm.chprobability[ findin(ChoicesAvailable(elm), action)[1] ]     # Record the prob with which the action was taken.
        newchoice = LevelFunction(elm, action)                                                 # What is the new level?
        elm.chprobability = HospUpdate(elm, newchoice)                                         # What are the new probabilities, given the new level?
        elm.level = newchoice                                                                  # Set the level to be the new choice.
        elm.levelhistory[i] = newchoice
        if el.fipscode == 48453
          println("    ", elm.fid, " ", elm.level, " ", round(elm.wtphist.w385[i],2))
          println("    ", patients.zips[78702].pdetutils)
        end 
      end
      # It would make sense to call the WTP update right here, in case it is not called earlier.
      
    end
    UpdateDeterministic(pats)                                                                  # Updates deterministic component of utility
    CleanWTPDict(wtpd1)                                                                        # cleans the WTP Dictionary
  end
  # TODO: Why return this?  Why not modify in place? this would require further modifications below 
  return Tex                                                                                   # Returns the whole state so the results can be written out.
end


"""
`EqAction`
NOT COMPLETE. NOT EXPORTED.
"""
function EqAction{T<:Fac}( H::T,i::Int64)
  action = StatsBase.sample( ChoicesAvailable(H), H.chprobability  )                 # Take the action
  # TODO - surely this next line can be fixed.  
  H.probhistory[i] = H.chprobability[ findin(ChoicesAvailable(H), action)[1] ]     # Record the prob with which the action was taken.
  newchoice = LevelFunction(H, action)                                                 # What is the new level?
  H.chprobability = HospUpdate(H, newchoice)                                         # What are the new probabilities, given the new level?
  H.level = newchoice                                                                  # Set the level to be the new choice.
  H.levelhistory[i] = newchoice
end





"""
`EntryProcess(Mkt::Market, i::Int64, T::Int64)`
Puts the entry process in a separate function.  
- for each mkt (county), draw an entrant: no entry, lev 1, 2 or 3
- Probs as given.  
- i is the current period in the sim - want to know where we are.
- T is the total length of each sim - want to know how long we might have to go.  

"""
function EntryProcess(el::Market, i::Int64, T::Int64)
  const entrants::Array{Int64,1} = [0, 1, 2, 3] 
  const entryprobs::Array{Float64,1} = [0.9895, 0.008, 0.0005, 0.002]
  entrant = StatsBase.sample(entrants, StatsBase.Weights(entryprobs))
  if entrant != 0
    entloc = NewEntrantLocation(el)                                                        # called on the market
    newfid = -floor(rand()*1e6)-1000000                                                    # all entrant fids negative to facilitate their removal later.
    entr = hospital( newfid, entloc[1], entloc[2], " Entrant $newfid ", el.fipscode, entrant, initial(entrant), Array{Int64,1}(T),
                     DemandHistory( Array{Int64,1}(T),  Array{Int64,1}(T), Array{Int64,1}(T), Array{Int64,1}(T), Array{Int64,1}(T),  Array{Int64,1}(T), Array{Int64,1}(T) ),
                     DemandHistory( Array{Int64,1}(T),  Array{Int64,1}(T), Array{Int64,1}(T), Array{Int64,1}(T), Array{Int64,1}(T),  Array{Int64,1}(T), Array{Int64,1}(T) ),
                     WTP( Array{Float64,1}(T),  Array{Float64,1}(T), Array{Float64,1}(T), Array{Float64,1}(T), Array{Float64,1}(T),  Array{Float64,1}(T), Array{Float64,1}(T) ),
                     StatsBase.Weights([0.1, 0.1, 0.1, 0.1]), Array{Float64,1}(T), neighbors(0, 0, 0, 0, 0, 0, 0, 0, 0), Array{Int64, 1}(), 0, false)
    entr.levelhistory[i] = entrant
    push!(el.config, entr)                                                                 # need to create a new record for this hospital in the market
    el.collection[newfid] = entr
    for elm in el.config                                                                   # need to add it to the dictionary too:
      NeighborAppend(elm, entr)
      NeighborAppend(entr, elm)
    end
    HospUpdate(entr, entrant)                                                             #"entrant" is the level
  end
end 






"""
`Termination(EmTex::EntireState)`
Takes an entire state (or the empty state for data recording) and returns "true" when every facility has been perturbed.
"""
function Termination(EmTex::EntireState)
  isdone = true
  for mrk in keys(EmTex.mkts) # iterates over markets
    # TODO - remove this comprehension - better to just iterate with & over all elements rather than allocating this.  
    isdone = (isdone)&(reduce(&, [ EmTex.mkts[mrk].noneqrecord[i] for i in keys(EmTex.mkts[mrk].noneqrecord) ] ))
  end
  return isdone
end

"""
`function RecordCopy{T<:ProjectModule.Fac}(ES::EntireState, h::T)`
This function will copy all the elements of the perturbed hospital
to the output record.
"""
function RecordCopy{T<:ProjectModule.Fac}(ES::EntireState, h::T)
  fips = ES.fipsdirectory[h.fid]
  for i = 1:length(h.probhistory)
    # probhistory
    ES.mkts[fips].collection[h.fid].probhistory[i] = h.probhistory[i]
    # WTP history
    ES.mkts[fips].collection[h.fid].wtphist.w385[i] = h.wtphist.w385[i]
    ES.mkts[fips].collection[h.fid].wtphist.w386[i] = h.wtphist.w386[i]
    ES.mkts[fips].collection[h.fid].wtphist.w387[i] = h.wtphist.w387[i]
    ES.mkts[fips].collection[h.fid].wtphist.w388[i] = h.wtphist.w388[i]
    ES.mkts[fips].collection[h.fid].wtphist.w389[i] = h.wtphist.w389[i]
    ES.mkts[fips].collection[h.fid].wtphist.w390[i] = h.wtphist.w390[i]
    ES.mkts[fips].collection[h.fid].wtphist.w391[i] = h.wtphist.w391[i]
    # levelhistory
    ES.mkts[fips].collection[h.fid].levelhistory[i] = h.levelhistory[i]
    # pdemandhist
    ES.mkts[fips].collection[h.fid].pdemandhist.demand385[i] = h.pdemandhist.demand385[i]
    ES.mkts[fips].collection[h.fid].pdemandhist.demand386[i] = h.pdemandhist.demand386[i]
    ES.mkts[fips].collection[h.fid].pdemandhist.demand387[i] = h.pdemandhist.demand387[i]
    ES.mkts[fips].collection[h.fid].pdemandhist.demand388[i] = h.pdemandhist.demand388[i]
    ES.mkts[fips].collection[h.fid].pdemandhist.demand389[i] = h.pdemandhist.demand389[i]
    ES.mkts[fips].collection[h.fid].pdemandhist.demand390[i] = h.pdemandhist.demand390[i]
    ES.mkts[fips].collection[h.fid].pdemandhist.demand391[i] = h.pdemandhist.demand391[i]
    #mdemandhist
    ES.mkts[fips].collection[h.fid].mdemandhist.demand385[i] = h.mdemandhist.demand385[i]
    ES.mkts[fips].collection[h.fid].mdemandhist.demand386[i] = h.mdemandhist.demand386[i]
    ES.mkts[fips].collection[h.fid].mdemandhist.demand387[i] = h.mdemandhist.demand387[i]
    ES.mkts[fips].collection[h.fid].mdemandhist.demand388[i] = h.mdemandhist.demand388[i]
    ES.mkts[fips].collection[h.fid].mdemandhist.demand389[i] = h.mdemandhist.demand389[i]
    ES.mkts[fips].collection[h.fid].mdemandhist.demand390[i] = h.mdemandhist.demand390[i]
    ES.mkts[fips].collection[h.fid].mdemandhist.demand391[i] = h.mdemandhist.demand391[i]
  end
end



"""
`PSim(T::Int64 ; di = ProjectModule.alldists, fi = fips, entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002])`
Runs a perturbed simulation - for each market, while there are hospitals I have not perturbed, runs a sim with one perturbed and the rest not.
The results are stored in EmptyState, which is an EntireState record instance.
This is for sure the slowest thing around.  
# for various benchmarking tasks.  
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
NewSim(50, Texas, patients);

Stampede:
@time PSim(1);  # 89.929505 seconds (65.28 M allocations: 1.263 GiB, 1.77% gc time)   x5.56
@time PSim(10); # 486.703537 seconds (314.47 M allocations: 5.919 GiB, 1.11% gc time) x4.86
@time PSim(40); # 1766.959950 seconds (1.16 G allocations: 21.890 GiB, 1.18% gc time) x4.38 (about 0.5 hours)
@time CombinedSim(1; T1 = 40);

Personal:
@time PSim(1);                   # 16.223982 seconds (65.29 M allocations: 1.263 GiB, 3.50% gc time)
@time PSim(10);                  # 100.210398 seconds (314.59 M allocations: 5.920 GiB, 1.96% gc time)
@time PSim(40);                  # 403.344224 seconds (1.16 G allocations: 21.912 GiB, 1.85% gc time)  
@time CombinedSim(1; T1 = 40);   # 415.258146 seconds (1.20 G allocations: 22.750 GiB, 1.90% gc time)

"""
function PSim(T::Int64; di = ProjectModule.alldists, fi = ProjectModule.fips)                           # fi = fips,
  const entrants::Array{Int64,1} = [0, 1, 2, 3]
  const entryprobs::Array{Float64,1} = [0.9895, 0.008, 0.0005, 0.002]
  EmptyState = CreateEmpty(fi, di, T);                                                                  # This is just a container of EntireState type - does not need linking.
  termflag = true                                                                                       # Initializes the termination flag.
  counter = 1
  arry1 = zeros(Int64, 1550)                                                                            # allocates an array for use in GenP.  Can be re-used.
  arry2 = zeros(Int64, 1550)                                                                            # allocates an array for use in GenM.  Can be re-used.
  wtparr1 = zeros(2,12)                                                                        # temporary array for CalcWTP
  Tex = CreateEmpty(fi, di, T)                                                                          # NB: New state once.
  pats = NewPatients(Tex);                                                                              # NB: New patient collection, linked to the new state.  Must be created AFTER "Tex."
  wtpd1 = WTPDict(Tex)                                                                         # creates a dict for WTPMap
  while termflag                                                                                        # true if there is some hospital which has not been perturbed.
    currentfac = Dict{Int64, Int64}()                                                                   # Dict{FID, fipscode} = {key, value}
    Restore(Tex)                                                                                        # cleans out all recorded data.
    UpdateDeterministic(pats)                                                                           # restores patient data to original state by updating utilities
    d1 = NewHospDict(Tex)                                                                               # creates a dict for GenP below.
    d2 = NewHospDict(Tex)                                                                               # creates a dict for GenM below
    for el in keys(EmptyState.mkts)
      if !reduce(&, [ EmptyState.mkts[el].noneqrecord[i] for i in keys(EmptyState.mkts[el].noneqrecord)]) # this is pretty fast: 0.000005 seconds (10 allocations: 240 bytes)
        # TODO - get rid of this comprehension.
        pfids = prod(hcat( [ [i, !EmptyState.mkts[el].noneqrecord[i]] for i in keys(EmptyState.mkts[el].noneqrecord) ]...) , 1) # 0.033497 seconds (14.66 k allocations: 794.414 KiB)
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
      # TODO - check that WTP gets updated in this function, just in case.  
      WTPMap(pats, Tex, wtpd1, wtparr1)
      WriteWTP(wtpd1, Tex, i)
      GenPChoices(pats, d1, arry1)                                                                     # this now modifies the dictionary in-place
      PDemandMap(d1, Tex, i)                                                                           # and this now cleans the dictionary up at the end, setting all demands to 0.
      GenMChoices(pats, d2, arry2)                                                                     # this now modifies the dictionary in-place
      MDemandMap(d2, Tex, i)                                                                           # and this now cleans the dictionary up at the end, setting all demands to 0.
      for el in Tex.ms
        if in(el.fipscode, pmarkets)                                                                   #NB: in( collection, element) !!
          EntryProcess(el, i, T)                                                                       # calls the EntryProcess - i.e., adds an entrant or not.  
          for elm in el.config
             if !elm.perturbed                                                                         # not perturbed, i.e., "perturbed" == false
               action = StatsBase.sample( ChoicesAvailable(elm), elm.chprobability )                   # Take the action
               # TODO - here add update of WTP given an action not equal to 10.
               # FIXME 
               elm.probhistory[i]= elm.chprobability[ findin(ChoicesAvailable(elm), action)[1] ]       # Record the prob with which the action was taken.
               newchoice = LevelFunction(elm, action)                                                  # What is the new level?
               elm.chprobability = HospUpdate(elm, newchoice)                                          # What are the new probabilities, given the new level?
               elm.level = newchoice                                                                   # Set the level to be the new choice.
               elm.levelhistory[i]=newchoice
             else # perturbed.
               action = StatsBase.sample( ChoicesAvailable(elm), HospPerturb(elm, elm.level,0.05))
               # TODO - here add update of WTP given an action not equal to 10.
               # FIXME 
               elm.probhistory[i]=elm.chprobability[findin(ChoicesAvailable(elm), action)[1]]
               newchoice = LevelFunction(elm, action)
               elm.chprobability = HospUpdate(elm, newchoice)
               elm.level = newchoice
               elm.levelhistory[i]=newchoice
             end
            # It would make sense to call the WTP update right here, in case it is not called earlier.
          end
        end
      end
      CleanWTPDict(wtpd1)                                                                        # cleans the WTP Dictionary
      UpdateDeterministic(pats)                                                                       # Updates deterministic component of utility for all patients and zips.
    end
    for fips in pmarkets                                                                              # a collection of fips codes
      RecordCopy(EmptyState, Tex.mkts[fips].collection[currentfac[fips]])                             # copies the record 
      EmptyState.mkts[fips].noneqrecord[ Tex.mkts[fips].collection[currentfac[fips]].fid ] = true     # update the value in the non-equilibrium record sim to true.
        #The next piece is not strictly necessary.  The config and collection point to the same underlying objects.  
      for nmb in 1:size(EmptyState.mkts[fips].config,1)                                               # iterate over the market config, which is an array.
        if EmptyState.mkts[fips].config[nmb].fid == Tex.mkts[fips].collection[currentfac[fips]].fid   # check for equality in the fids
          RecordCopy(EmptyState, Tex.mkts[fips].collection[currentfac[fips]])                         # here I am reassigning.  
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

This produces a series of matrices: 4 of them.  Three of them have three rows and #DRGs columns.  
- private is a count of the number of privately insured patients admitted at each of seven DRGS
- medicaid is the same for the medicaid patients - a count at seven DRGs 
- wtp_out computes WTP measures by DRG.
- the rows correspond to different levels, so one row is the history of demand for seven drg's while the level was 2.

- the fourth thing returned is transitions, which is a record of transitions across levels.

Prior to output, the shape of private, medicaid and WTP is: 
[385_1, 385_2, 385_3;
 386_1, 386_2, 386_3;
 387_1, 387_2, 387_3;
 388_1, 388_2, 388_3;
 389_1, 389_2, 389_3;
 390_1, 390_2, 390_3;
 391_1, 391_2, 391_3; ]

where these entries are demand/WTP_level.

Then at the end this is reshaped so that the form is: [385_1 ... 391_1 385_2 ... 391_2 385_3 ... 391_3]

Transitions is: 9x1. For the form of that, see TransitionGen above.
[ ; 

  ]

Output is a 4-tuple of vectors: private, medicaid, wtp, transitions.
This is given to ResultsOut() below.  

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
- Takes the output of CondSum, which is a 4-tuple of vectors private, medicaid, wtp_out, transitions
- Assigns these elements to the 1x40 array arr.
- Multiplies by the discount rate.  FIXME - this is not quite right now. 
- Writes out the results as a large matrix, with one row for each firm.  


  Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
  Nexas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
  Restore(Texas) # assign all to 0 or 1
  Restore(Nexas)

  # TODO - this should definitely have a test.  

"""
function ResultsOut(Tex::EntireState, OtherTex::EntireState; T::Int64 = 50, beta::Float64 = 0.95,  dim2::Int64 = 81) #dim2 - 33 paramsx2 + 7x2 records of medicaid volumes + one identifying FID
  const drgamt::Array{Float64,1} = [12038.83, 66143.19, 19799.52, 4044.67, 6242.39, 1329.98, 412.04]
  dim1 = Tex.fipsdirectory.count
  outp = Array{Float64,2}(dim1, dim2)
  fids = [k for k in keys(Tex.fipsdirectory)]
  disc::Float64 = (1-beta^(T+1))/(1-beta) 
  for el in 1:size(fids,1)
    outp[el,1] = fids[el]                                                           # Write out all of the fids as an ID in the first column.
  end
  for el in keys(Tex.fipsdirectory)                                                 # Now this is all of the hospitals by FID
    hosp = Tex.mkts[Tex.fipsdirectory[el]].collection[el]
    outprob = prod(hosp.probhistory)                                                # Prob of the outcome.
    private, medicaid, wtp_out, transitions = CondSum(hosp)
    arr = zeros(1, 40)
    arr[1] = disc*outprob*dot(wtp_out[1:7], private[1:7])                       # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 1
    arr[2] = disc*outprob*dot(wtp_out[8:14], private[8:14])                     # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 2
    arr[3] = disc*outprob*dot(wtp_out[15:21], private[15:21])                   # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 3
    # The next lines are patients summed over types.  Costs are treated as the same over Medicaid and privately insured.
    arr[4] = disc*outprob*(private[1]+medicaid[1])                              # 385_p_1 + 385_m_1                             
    arr[5] = disc*outprob*(private[8]+medicaid[8])                              # 385_p_2 + 385_m_2
    arr[6] = disc*outprob*(private[15]+medicaid[15])                            # 385_p_3 + 385_m_3 
    arr[7] = disc*outprob*(private[2]+medicaid[2])                              # 386_p_1 + 386_m_1 
    arr[8] = disc*outprob*(private[9]+medicaid[9])                              # 386_p_2 + 386_m_2
    arr[9] = disc*outprob*(private[16]+medicaid[16])                            # 386_p_3 + 386_m_3
    arr[10] = disc*outprob*(private[3]+medicaid[3])                             # 387_p_1 + 387_m_1
    arr[11] = disc*outprob*(private[10]+medicaid[10])                           # 387_p_2 + 387_m_2
    arr[12] = disc*outprob*(private[17]+medicaid[17])                           # 387_p_3 + 387_m_3
    arr[13] = disc*outprob*(private[4]+medicaid[4])                             # 388_p_1 + 388_m_1
    arr[14] = disc*outprob*(private[11]+medicaid[11])                           # 388_p_2 + 388_m_2
    arr[15] = disc*outprob*(private[18]+medicaid[18])                           # 388_p_3 + 388_m_3
    arr[16] = disc*outprob*(private[5]+medicaid[5])                             # 389_p_1 + 389_m_1
    arr[17] = disc*outprob*(private[12]+medicaid[12])                           # 389_p_2 + 389_m_2
    arr[18] = disc*outprob*(private[19]+medicaid[19])                           # 389_p_3 + 389_m_3
    arr[19] = disc*outprob*(private[6]+medicaid[6])                             # 390_p_1 + 390_m_1
    arr[20] = disc*outprob*(private[13]+medicaid[13])                           # 390_p_2 + 390_m_2
    arr[21] = disc*outprob*(private[20]+medicaid[20])                           # 390_p_3 + 390_m_3
    arr[22] = disc*outprob*(private[7]+medicaid[7])                             # 391_p_1 + 391_m_1
    arr[23] = disc*outprob*(private[14]+medicaid[14])                           # 391_p_2 + 391_m_2
    arr[24] = disc*outprob*(private[21]+medicaid[21])                           # 391_p_3 + 391_m_3
    arr[25] = disc*outprob*(medicaid[1]+medicaid[8]+medicaid[15])*drgamt[1]     # (385_m_1 + 385_m_2 + 385_m_3)* revenue avg. at DRG 385
    arr[26] = disc*outprob*(medicaid[2]+medicaid[9]+medicaid[16])*drgamt[2]     # (386_m_1 + 386_m_2 + 386_m_3)* revenue avg. at DRG 386
    arr[27] = disc*outprob*(medicaid[3]+medicaid[10]+medicaid[17])*drgamt[3]    # (387_m_1 + 387_m_2 + 387_m_3)* revenue avg. at DRG 387
    arr[28] = disc*outprob*(medicaid[4]+medicaid[11]+medicaid[18])*drgamt[4]    # (388_m_1 + 388_m_2 + 388_m_3)* revenue avg. at DRG 388
    arr[29] = disc*outprob*(medicaid[5]+medicaid[12]+medicaid[19])*drgamt[5]    # (389_m_1 + 389_m_2 + 389_m_3)* revenue avg. at DRG 389
    arr[30] = disc*outprob*(medicaid[6]+medicaid[13]+medicaid[20])*drgamt[6]    # (390_m_1 + 390_m_2 + 390_m_3)* revenue avg. at DRG 390
    arr[31] = disc*outprob*(medicaid[7]+medicaid[14]+medicaid[21])*drgamt[7]    # (391_m_1 + 391_m_2 + 391_m_3)* revenue avg. at DRG 391
    arr[32:end] = (alltrans = disc*outprob*transitions)                         # Transitions.
    index = findfirst(outp[:,1], hosp.fid)                                          # find where the fid is in the list.
    outp[index, 2:41] = arr
    # NB: Here starts the second state record.
    hosp_neq = OtherTex.mkts[OtherTex.fipsdirectory[el]].collection[el]             # Find the record in the OTHER EntireState
    outprobn = prod(hosp_neq.probhistory)                                           # Prob of the outcome.
    privaten, medicaidn, wtp_outn, transitionsn = CondSum(hosp_neq)
    narr = zeros(1, 40)
    narr[1] = disc*outprobn*dot(wtp_outn[1:7], privaten[1:7])                     # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 1
    narr[2] = disc*outprobn*dot(wtp_outn[8:14], privaten[8:14])                   # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 2
    narr[3] = disc*outprobn*dot(wtp_outn[15:21], privaten[15:21])                 # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 3
    # The next lines are patients summed over types.  Costs are treated as the same over Medicaid and privately insured.
    narr[4] = disc*outprobn*(privaten[1]+medicaidn[1])                            # 385_p_1 + 385_m_1                   
    narr[5] = disc*outprobn*(privaten[8]+medicaidn[8])                            # 385_p_2 + 385_m_2
    narr[6] = disc*outprobn*(privaten[15]+medicaidn[15])                          # 385_p_3 + 385_m_3
    narr[7] = disc*outprobn*(privaten[2]+medicaidn[2])                            # 386_p_1 + 386_m_1
    narr[8] = disc*outprobn*(privaten[9]+medicaidn[9])                            # 386_p_2 + 386_m_2
    narr[9] = disc*outprobn*(privaten[16]+medicaidn[16])                          # 386_p_3 + 386_m_3
    narr[10] = disc*outprobn*(privaten[3]+medicaidn[3])                           # 387_p_1 + 387_m_1
    narr[11] = disc*outprobn*(privaten[10]+medicaidn[10])                         # 387_p_2 + 387_m_2
    narr[12] = disc*outprobn*(privaten[17]+medicaidn[17])                         # 387_p_3 + 387_m_3
    narr[13] = disc*outprobn*(privaten[4]+medicaidn[4])                           # 388_p_1 + 388_m_1
    narr[14] = disc*outprobn*(privaten[11]+medicaidn[11])                         # 388_p_2 + 388_m_2
    narr[15] = disc*outprobn*(privaten[18]+medicaidn[18])                         # 388_p_3 + 388_m_3
    narr[16] = disc*outprobn*(privaten[5]+medicaidn[5])                           # 389_p_1 + 389_m_1
    narr[17] = disc*outprobn*(privaten[12]+medicaidn[12])                         # 389_p_2 + 389_m_2
    narr[18] = disc*outprobn*(privaten[19]+medicaidn[19])                         # 389_p_3 + 389_m_3
    narr[19] = disc*outprobn*(privaten[6]+medicaidn[6])                           # 390_p_1 + 390_m_1
    narr[20] = disc*outprobn*(privaten[13]+medicaidn[13])                         # 390_p_2 + 390_m_2
    narr[21] = disc*outprobn*(privaten[20]+medicaidn[20])                         # 390_p_3 + 390_m_3
    narr[22] = disc*outprobn*(privaten[7]+medicaidn[7])                           # 391_p_1 + 391_m_1
    narr[23] = disc*outprobn*(privaten[14]+medicaidn[14])                         # 391_p_2 + 391_m_2
    narr[24] = disc*outprobn*(privaten[21]+medicaidn[21])                         # 391_p_3 + 391_m_3
    narr[25] = disc*outprobn*(medicaidn[1]+medicaidn[8]+medicaidn[15])*drgamt[1]  # (385_m_1 + 385_m_2 + 385_m_3)* revenue avg. at DRG 385
    narr[26] = disc*outprobn*(medicaidn[2]+medicaidn[9]+medicaidn[16])*drgamt[2]  # (386_m_1 + 386_m_2 + 386_m_3)* revenue avg. at DRG 386
    narr[27] = disc*outprobn*(medicaidn[3]+medicaidn[10]+medicaidn[17])*drgamt[3] # (387_m_1 + 387_m_2 + 387_m_3)* revenue avg. at DRG 387
    narr[28] = disc*outprobn*(medicaidn[4]+medicaidn[11]+medicaidn[18])*drgamt[4] # (388_m_1 + 388_m_2 + 388_m_3)* revenue avg. at DRG 388
    narr[29] = disc*outprobn*(medicaidn[5]+medicaidn[12]+medicaidn[19])*drgamt[5] # (389_m_1 + 389_m_2 + 389_m_3)* revenue avg. at DRG 389
    narr[30] = disc*outprobn*(medicaidn[6]+medicaidn[13]+medicaidn[20])*drgamt[6] # (390_m_1 + 390_m_2 + 390_m_3)* revenue avg. at DRG 390
    narr[31] = disc*outprobn*(medicaidn[7]+medicaidn[14]+medicaidn[21])*drgamt[7] # (391_m_1 + 391_m_2 + 391_m_3)* revenue avg. at DRG 391
    narr[32:end] = disc*outprobn*transitionsn                                 # Transitions - 9 of them.
    outp[index, 42:end] = narr
  end
  return sortrows(outp, by=x->x[1])                                                    # sort by first column (fid)
end




"""
`ResultsOutVariant`

This function will compute the results of the simulation on the assumption that there 
are cost parameters which vary only by level, not by DRG and level.  

So instead of assuming that the profit is...

∑ βᵗ α WTP() - γ₋{drg, level} Vᵖ  

the profit is 

∑ βᵗ α WTP() - γ₋{level} Vᵖ



"""
function ResultsOutVariant(Tex::EntireState, OtherTex::EntireState; T::Int64 = 50, beta::Float64 = 0.95) 
  # there are ultimately fewer parameters by about... 6?  or 12?  
  # how many exactly... WTP1, WTP2, WTP3, COST1, COST2, COST3, REV385, REV386, REV387, REV388, REV389, REV390, REV391, Transitionsx9
  # this should be 22
  # TODO - also need to add some kind of composite for the mothers...?  Or what?  
  const params::Int64 = 22
  const dim2::Int64 = 45
  const drgamt::Array{Float64,1} = [12038.83, 66143.19, 19799.52, 4044.67, 6242.39, 1329.98, 412.04]
  const weight385::Float64 = 1.38
  const weight386::Float64 = 4.57
  const weight387::Float64 = 3.12
  const weight388::Float64 = 1.88
  const weight389::Float64 = 3.20
  const weight390::Float64 = 1.13
  const weight391::Float64 = 0.15
  dim1 = Tex.fipsdirectory.count
  outp = Array{Float64,2}(dim1, dim2)
  fids = [k for k in keys(Tex.fipsdirectory)]
  disc::Float64 = (1-beta^(T+1))/(1-beta)
  for el in 1:size(fids,1)
    outp[el,1] = fids[el]                                                           # Write out all of the fids as an ID in the first column.
  end
  for el in keys(Tex.fipsdirectory)                                                 # Now this is all of the hospitals by FID
    hosp = Tex.mkts[Tex.fipsdirectory[el]].collection[el]
    outprob = prod(hosp.probhistory)                                                # Prob of the outcome.
    private, medicaid, wtp_out, transitions = CondSum(hosp)
    arr = zeros(1, params)
    arr[1] = disc*outprob*dot(wtp_out[1:7], private[1:7])                       # This is WTP over all DRGS * private patient vols at corresponding DRG, over all periods at level 1
    arr[2] = disc*outprob*dot(wtp_out[8:14], private[8:14])                     # This is WTP over all DRGS * private patient vols at corresponding DRG, over all periods at level 2
    arr[3] = disc*outprob*dot(wtp_out[15:21], private[15:21])                   # This is WTP over all DRGS * private patient vols at corresponding DRG, over all periods at level 3
    # this is the combined cost at level 1, multiplied by weights: (385_m_1 + 385_p_1)*weight385 + ... + (391_m_1+391_p_1)*weight391 
    arr[4] = disc*outprob*(weight385*(private[1]+medicaid[1])+weight386*(private[2]+medicaid[2])+weight387*(private[3]+medicaid[3])+weight388*(private[4]+medicaid[4])+weight389*(private[5]+medicaid[5])+weight390*(private[6]+medicaid[6])+weight391*(private[7]+medicaid[7])) 
    # this is the combined cost at level 2, multiplied by weights: (385_m_2 + 385_p_2)*weight385 + ... + (391_m_2+391_p_2)*weight391 
    arr[5] = disc*outprob*(weight385*(private[8]+medicaid[8])+weight386*(private[9]+medicaid[9])+weight387*(private[10]+medicaid[10])+weight388*(private[11]+medicaid[11])+weight389*(private[12]+medicaid[12])+weight390*(private[13]+medicaid[13])+weight391*(private[14]+medicaid[14]))
    # this is the combined cost at level 3, multiplied by weights: (385_m_3 + 385_p_3)*weight385 + ... + (391_m_3+391_p_3)*weight391 
    arr[6] = disc*outprob*(weight385*(private[15]+medicaid[15])+weight386*(private[16]+medicaid[16])+weight387*(private[17]+medicaid[17])+weight388*(private[18]+medicaid[18])+weight389*(private[19]+medicaid[19])+weight390*(private[20]+medicaid[20])+weight391*(private[21]+medicaid[21]))
    arr[7] = disc*outprob*(medicaid[1]+medicaid[8]+medicaid[15])*drgamt[1]    # Patients*revenue avg. at DRG 385
    arr[8] = disc*outprob*(medicaid[2]+medicaid[9]+medicaid[16])*drgamt[2]    # Patients*revenue avg. at DRG 386
    arr[9] = disc*outprob*(medicaid[3]+medicaid[10]+medicaid[17])*drgamt[3]   # Patients*revenue avg. at DRG 387
    arr[10] = disc*outprob*(medicaid[4]+medicaid[11]+medicaid[18])*drgamt[4]  # Patients*revenue avg. at DRG 388
    arr[11] = disc*outprob*(medicaid[5]+medicaid[12]+medicaid[19])*drgamt[5]  # Patients*revenue avg. at DRG 389
    arr[12] = disc*outprob*(medicaid[6]+medicaid[13]+medicaid[20])*drgamt[6]  # Patients*revenue avg. at DRG 390
    arr[13] = disc*outprob*(medicaid[7]+medicaid[14]+medicaid[21])*drgamt[7]  # Patients*revenue avg. at DRG 391, 13 parameters to here.
    arr[14:end] = (alltrans = (beta^T)*outprob*transitions)                   # 9 here. 
    index = findfirst(outp[:,1], hosp.fid)                                    # find where the fid is in the list.
    outp[index, 2:(params+1)] = arr
    # NB: Here starts the second state record.
    hosp_neq = OtherTex.mkts[OtherTex.fipsdirectory[el]].collection[el]       # Find the record in the OTHER EntireState
    outprobn = prod(hosp_neq.probhistory)                                     # Prob of the outcome.
    privaten, medicaidn, wtp_outn, transitionsn = CondSum(hosp_neq)
    narr = zeros(1, 22)
    narr[1] = disc*outprobn*dot(wtp_outn[1:7], privaten[1:7])             # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 1
    narr[2] = disc*outprobn*dot(wtp_outn[8:14], privaten[8:14])           # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 2
    narr[3] = disc*outprobn*dot(wtp_outn[15:21], privaten[15:21])         # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 3
    narr[4] = disc*outprobn*(weight385*(privaten[1]+medicaidn[1])+weight386*(privaten[2]+medicaidn[2])+weight387*(privaten[3]+medicaidn[3])+weight388*(privaten[4]+medicaidn[4])+weight389*(privaten[5]+medicaidn[5])+weight390*(privaten[6]+medicaidn[6])+weight391*(privaten[7]+medicaidn[7]))                   
    narr[5] = disc*outprobn*(weight385*(privaten[8]+medicaidn[8])+weight386*(privaten[9]+medicaidn[9])+weight387*(privaten[10]+medicaidn[10])+weight388*(privaten[11]+medicaidn[11])+weight389*(privaten[12]+medicaidn[12])+weight390*(privaten[13]+medicaidn[13])+weight391*(privaten[14]+medicaidn[14]))
    narr[6] = disc*outprobn*(weight385*(privaten[15]+medicaidn[15])+weight386*(privaten[16]+medicaidn[16])+weight387*(privaten[17]+medicaidn[17])+weight388*(privaten[18]+medicaidn[18])+weight389*(privaten[19]+medicaidn[19])+weight390*(privaten[20]+medicaidn[20])+weight391*(privaten[21]+medicaidn[21]))
    narr[7] = disc*outprobn*(medicaidn[1]+medicaidn[8]+medicaidn[15])*drgamt[1]
    narr[8] = disc*outprobn*(medicaidn[2]+medicaidn[9]+medicaidn[16])*drgamt[2]
    narr[9] = disc*outprobn*(medicaidn[3]+medicaidn[10]+medicaidn[17])*drgamt[3]
    narr[10] = disc*outprobn*(medicaidn[4]+medicaidn[11]+medicaidn[18])*drgamt[4]
    narr[11] = disc*outprobn*(medicaidn[5]+medicaidn[12]+medicaidn[19])*drgamt[5]
    narr[12] = disc*outprobn*(medicaidn[6]+medicaidn[13]+medicaidn[20])*drgamt[6]
    narr[13] = disc*outprobn*(medicaidn[7]+medicaidn[14]+medicaidn[21])*drgamt[7]
    narr[14:end] = disc*outprobn*transitionsn
    outp[index, (params+2):dim2] = narr
  end
  return sortrows(outp, by=x->x[1])                                                    # sort by first column (fid)
end 





"""
`OuterSim(MCcount::Int; T1::Int64 = 3, fi = ProjectModule.fips, di=ProjectModule.alldists)`
Runs the Monte Carlo - Equilibrium and Non-equilibrium simulations for each market MCcount times.
Note that the reduction is (+), but that includes adding the fids, so this must be divided by MCcount
to return correct results.

Testing:
OuterSim(3)

# TODO - need to fix this to return a tuple.  that should be possible. 
# TODO - how much time per loop does remaking patients and CreateEmpty take???  about 0.17 seconds.  Not too much.  

"""
function OuterSim(MCcount::Int; T1::Int64 = 3, fi = ProjectModule.fips, di = ProjectModule.alldists)
  outp = @sync @parallel (+) for j = 1:MCcount
    println("Current iteration ", j)
    TexasEq = CreateEmpty(fi, di, T1)
    #TexasNeq = MakeNew(fi, da);                                                                         # Returns a separate EntireState.
    eq_patients = NewPatients(TexasEq)                                                                   # Separate patients - these linked to Eq Entire State.
    #neq_patients = NewPatients(TexasNeq)                                                                # Separate patients - these linked to Neq Entire State.
    ResultsOut(NewSim(T1, TexasEq, eq_patients), PSim(T1); T = T1)                                       # simulates and writes out the results.
  end
  outp[:,1] = outp[:,1]/MCcount                                                                          # Combined by (+) so reproduce the fids by dividing.
  return outp
end



"""
`CombinedSim(MCcount::Int; T1::Int64 = 3, fi = ProjectModule.fips, di = ProjectModule.alldists)` 
Goal here to get two outputs from one simulation.  

TODO - something is making this REALLY slow on stampede2.  REALLY slow.  What the hell?
NB - this comparison isn't quite fair yet because Stampede does not have the most recent version w/ changes to WTPMap etc.  
Stampede MKL: @time NewSim(5, Texas, patients); #6.970327 seconds (4.71 M allocations: 104.489 MiB, 1.93% gc time)
Stampede MKL: @time NewSim(40, Texas, patients); #55.357845 seconds (37.32 M allocations: 828.682 MiB, 1.37% gc time)

Personal: @time NewSim(5, Texas, patients); #1.751509 seconds (4.69 M allocations: 92.628 MiB, 2.81% gc time)
Personal: @time NewSim(40, Texas, patients); #13.861571 seconds (36.96 M allocations: 730.009 MiB, 1.82% gc time)


# STAMPEDE2 BENCHMARKS: 

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 5);
patients = NewPatients(Texas);
@time NewSim(5, Texas, patients); #9.312147 seconds


Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 40);
patients = NewPatients(Texas);
@time NewSim(40, Texas, patients); #73.896560 seconds (37.03 M allocations: 823.834 MiB, 1.02% gc time)
# if this takes 5x longer on stampede, then...
# PSim(20) = 210*5 = 1050 ≈ 17 mins...
# PSim(40) = 441*5 = 2205 seconds ≈ 0.5 hours... Why so long?  


PSim(1); # force compilation, though that's almost nothing in this case.  Or maybe not... ?


# COMPUTER Benchmarks 

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 20);
eq_patients = NewPatients(TexasEq);
@time NewSim(20, TexasEq, eq_patients); # 7 seconds.  
PSim(20) # - probably 30*7 = 210 seconds... roughly.  (Actually 214)


TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 40);
eq_patients = NewPatients(TexasEq);
@time NewSim(40, TexasEq, eq_patients); # 14 seconds @ 1 thread,  17 second @ 4 threads.  Weird.  

@time PSim(40); #  441.245888 seconds @ 1 thread, 466.358290 seconds @ 4 threads.  


What about DoubleResults?

DoubleResults(NewSim(1, TexasEq, eq_patients), PSim(1)); # compile 


TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 20);
eq_patients = NewPatients(TexasEq);
@time DoubleResults(NewSim(20, TexasEq, eq_patients), PSim(20)); # 223.502099 seconds - the sum of the individual times.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 40);
eq_patients = NewPatients(TexasEq);
@time DoubleResults(NewSim(40, TexasEq, eq_patients), PSim(40)) ; # 461.640865 seconds - basically the sum.


"""
function CombinedSim(MCcount::Int; T1::Int64 = 3, fi = ProjectModule.fips, di = ProjectModule.alldists)
  outp1, outp2 = @sync @parallel (ArrayTupleSum) for j = 1:MCcount
    println("iteration: ", j)
    TexasEq = CreateEmpty(fi, di, T1)
    eq_patients = NewPatients(TexasEq)
    # TODO - what does T do in the next function? 
    DoubleResults(NewSim(T1, TexasEq, eq_patients), PSim(T1); T = T1)
    # This will create a new TexasEq and patients every loop - can they be cleaned up instead?  
  end 
  outp1[:,1] = outp1[:,1]/MCcount 
  outp2[:,1] = outp2[:,1]/MCcount
  return outp1, outp2
end 



"""
`DoubleResults`
Simply calls ResultsOut and ResultsOutVariant on two state arguments, then returns the tuple of their outputs.
Designed to be called in CombinedSim.

### Testing ###

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
NewSim(50, Texas, patients);

DoubleResults(Texas, Texas)

"""
# TODO - what does T do here?  Nothing, it seems... 
function DoubleResults(Tex::EntireState, OtherTex::EntireState; T::Int64 = 50)
  # return two arrays, this can be combined with (+)
  # no - tuple of arrays cannot be combined because tuples are immutable.  
  return ResultsOut(Tex, OtherTex), ResultsOutVariant(Tex, OtherTex)
end 

"""
`ArrayTupleSum(t1::Tuple{Array{Float64,2},Array{Float64,2}}, t2::Tuple{Array{Float64,2},Array{Float64,2}})`

this is kind of a moronic function.  `DoubleResults` returns the tuple of arrays from `ResultsOut` and 
`ResultsOutVariant`.  Tuples are immutable, so this cannot be combined in a parallel-for in `CombinedSim`,
but it does work to define this stupid function.  It is not super efficient, but this part of the overhead is
probably small.  There are NO size checks or any other reasonable things to make sure it works, so only use 
it in that one place.  
"""
function ArrayTupleSum(t1::Tuple{Array{Float64,2},Array{Float64,2}}, t2::Tuple{Array{Float64,2},Array{Float64,2}})
  return t1[1]+t2[1], t1[2]+t2[2]
end 





## NB: To clean up


"""
`function CleanDemandHist(d::DemandHistory)`
Sets values back to zero.
"""
function CleanDemandHistory(d::DemandHistory)
  for i = 1:length(d.demand385)
    d.demand385[i] = 0
    d.demand386[i] = 0
    d.demand387[i] = 0
    d.demand388[i] = 0
    d.demand389[i] = 0
    d.demand390[i] = 0
    d.demand391[i] = 0
  end
end

"""
`function CleanLevelHistory(hos::hospital)`
Cleans out the level history AFTER the first one is written back.
"""
function CleanLevelHistory(hos::hospital)
  hos.levelhistory[1] = hos.init.level
  for i = 2:length(hos.levelhistory)
    hos.levelhistory[i] = 0
  end
end

"""
`function CleanWTPHistory`
Cleans out the WTP history.  Sets all to 0.0 since all are Float64.
"""
function CleanWTPHistory(hos::hospital)
  for i = 1:length(hos.wtphist.w385)
    hos.wtphist.w385[i] = 0.0
    hos.wtphist.w386[i] = 0.0
    hos.wtphist.w387[i] = 0.0
    hos.wtphist.w388[i] = 0.0
    hos.wtphist.w389[i] = 0.0
    hos.wtphist.w390[i] = 0.0
    hos.wtphist.w391[i] = 0.0
  end
end

"""
`function CleanProbHistory`
Sets prob history back to 0.0, except for the first entry which is 1.

TODO - do I want 0 or 1?  
"""
function CleanProbHistory(h::hospital)
  h.probhistory[1] = 1.0
  for i = 2:length(h.probhistory)
    # NB - changed to 1, to avoid 0 output.  FIXME.  check.  
    h.probhistory[i] = 1.0
  end
end

"""
`HospitalClean(hos::hospital)`
resets the hospital to the initial state after a run of the simulation.
some fields don't change: fid, lat, long, name, fips, bedcount,
eventually perturbed will probably change.
"""
function HospitalClean(hos::hospital)
  hos.level = hos.init.level                 #Reset level history to initial value
  CleanLevelHistory(hos)           #Set record length back to zero.
  CleanDemandHistory(hos.mdemandhist)  #empty demand history
  CleanDemandHistory(hos.pdemandhist)  #empty demand history
  CleanWTPHistory(hos)  #empty WTP
  #hos.chprobability = Weights([0]) # Reset this in Restore below when all states and levels have been reset.
  CleanProbHistory(hos)             #Empty history of choices.
  hos.neigh.level105 = 0
  hos.neigh.level205 = 0
  hos.neigh.level305 = 0
  hos.neigh.level1515 = 0
  hos.neigh.level2515 = 0
  hos.neigh.level3515 = 0
  hos.neigh.level11525 = 0
  hos.neigh.level21525 = 0
  hos.neigh.level31525 = 0
  hos.hood = Array{Int64,1}()
  hos.perturbed = false                            #For now this is always false.
end


"""
`RemoveEntrant(Tex::EntireState)`
Finds entrants and removes them.
Must go over ALL markets.
Remove from collection and config.
"""
function RemoveEntrant(Tex::EntireState)
  for k1 in keys(Tex.mkts)
    ents = Array{Int64,1}() # this allocates a little.
    for k2 in keys(Tex.mkts[k1].collection)
      if k2 < 0
        push!(ents, k2)
      end
    end
    if length(ents)>0 # nonzero length means more than one entrant.
      for el in ents
        delete!(Tex.mkts[k1].collection, el)
        deleteat!(Tex.mkts[k1].config, FidFindFirst(Tex.mkts[k1], el))
      end
    end
  end
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
  RemoveEntrant(Tex)
  NeighborFix(Tex) # Restores all neighbors to both hosp.neigh and hosp.hood.
  for mkt in Tex.ms
    for hos in mkt.config
      HospUpdate(hos, hos.level; update = true) # HospUpdate should now fix these.
    end
  end
end



###
