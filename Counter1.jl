# First Counterfactual.





# Create and fill state and Patient collection:
# Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
# CMakeIt(Tex, ProjectModule.fips);
# FillState(Tex, ProjectModule.alldists);
# patients = NewPatients(Tex);
#





#NB:  Remember to update the deterministic utilities once this part has been changed.




"""
`CategoryReminder(v::Int64)`
Reminds me what weight category the integers refer to
takes an integer and returns the weight interval.
"""
function CategoryReminder(v::Int64)
  if v ==1
    return "<= 499 grams"
  elseif  v ==2
    return "500-749 grams"
  elseif  v ==3
    return "750-999 grams"
  elseif  v ==4
    return "1000-1249 grams"
  elseif  v ==5
    return "1250-1499 grams"
  elseif  v ==6
    return "1500-1999 grams"
  elseif  v ==7
    return "2000-2499 grams"
  elseif  v ==8
    return "2500-2999 grams"
  elseif  v ==9
    return "3000-3499 grams"
  elseif  v ==10
    return "3500-3999 grams"
  elseif  v ==11
    return "4000-4499 grams"
  elseif  v ==12
    return "4500-4999 grams"
  elseif v ==13
    return "5000-8165 grams"
  else
    return "Not a valid option: 1-13 only"
  end
end






"""
`Payoff(ppats::Dict{Int64, ProjectModule.patientcount}, mpats::Dict{Int64, ProjectModule.patientcount}, Tex::EntireState, wtp::Dict{Int64,Float64}; params = [])`
Computes the actual firm payoffs.  Uses parameters computed from one run of the LTE.
"""
function Payoff(ppats::Dict{Int64, ProjectModule.patientcount}, mpats::Dict{Int64, ProjectModule.patientcount}, Tex::EntireState, wtp::Dict{Int64,Float64} ;
                alf1::Float64 = 29182.967,
                alf2::Float64 = 22167.6375,
                alf3::Float64 = 23074.8403,
                gamma_1_385::Float64 = 34628.8402,
                gamma_2_385::Float64 = 14921.003,
                gamma_3_385::Float64 = 12822.723,
                gamma_1_386::Float64 = 104578.867,
                gamma_2_386::Float64 = 95366.0004,
                gamma_3_386::Float64 = 69353.471,
                gamma_1_387::Float64 = 34498.5261,
                gamma_2_387::Float64 = 48900.8396,
                gamma_3_387::Float64 = 24639.0552,
                gamma_1_388::Float64 = 26561.8688,
                gamma_2_388::Float64 = 20895.5001,
                gamma_3_388::Float64 = 29775.8381,
                gamma_1_389::Float64 = 20653.5821, 
                gamma_2_389::Float64 = 20102.2097,
                gamma_3_389::Float64 = 8279.774,
                gamma_1_390::Float64 = 7372.3301,
                gamma_2_390::Float64 = 2514.8717,
                gamma_3_390::Float64 = 26113.4462,
                gamma_1_391::Float64 = 27018.9915,
                gamma_2_391::Float64 = 15079.2889,
                gamma_3_391::Float64 = 1912.7285 ) # params
  outp::Dict{Int64,Float64} = Dict{Int64,Float64}() # Dict(k => 0.0 for k in keys(ppats))
  for k1 in keys(ppats) # put the dictionary together.
    if k1 != 0 # don't add OO.
      outp[k1] = 0.0
    end 
  end 
  for k in keys(ppats) # multiplied by patient volumes or not?
    if k != 0 # don't include OO.
      if Tex.mkts[Tex.fipsdirectory[k]].collection[k].level == 1
        outp[k] = alf1*wtp[k]*(sum(ppats[k])+sum(mpats[k])) - gamma_1_385*ppats[k].count385 - gamma_1_385*mpats[k].count385 - gamma_1_386*ppats[k].count386 - gamma_1_386*mpats[k].count386 - gamma_1_387*ppats[k].count387 - gamma_1_387*mpats[k].count387 - gamma_1_388*mpats[k].count388 - gamma_1_388*ppats[k].count388 - gamma_1_389*mpats[k].count389 - gamma_1_389*ppats[k].count389 - gamma_1_390*ppats[k].count390 - gamma_1_390*mpats[k].count390 - gamma_1_391*ppats[k].count391 - gamma_1_391*mpats[k].count391
      elseif Tex.mkts[Tex.fipsdirectory[k]].collection[k].level == 2
        outp[k] = alf2*wtp[k]*(sum(ppats[k])+sum(mpats[k])) - gamma_2_385*ppats[k].count385 - gamma_2_385*mpats[k].count385 - gamma_2_386*ppats[k].count386 - gamma_2_386*mpats[k].count386 - gamma_2_387*ppats[k].count387 - gamma_2_387*mpats[k].count387 - gamma_2_388*mpats[k].count388 - gamma_2_388*ppats[k].count388 - gamma_2_389*mpats[k].count389 - gamma_2_389*ppats[k].count389 - gamma_2_390*ppats[k].count390 - gamma_2_390*mpats[k].count390 - gamma_2_391*ppats[k].count391 - gamma_2_391*mpats[k].count391
      else # level is 3
        outp[k] = alf3*wtp[k]*(sum(ppats[k])+sum(mpats[k])) - gamma_3_385*ppats[k].count385 - gamma_3_385*mpats[k].count385 - gamma_3_386*ppats[k].count386 - gamma_3_386*mpats[k].count386 - gamma_3_387*ppats[k].count387 - gamma_3_387*mpats[k].count387 - gamma_3_388*mpats[k].count388 - gamma_3_388*ppats[k].count388 - gamma_3_389*mpats[k].count389 - gamma_3_389*ppats[k].count389 - gamma_3_390*ppats[k].count390 - gamma_3_390*mpats[k].count390 - gamma_3_391*ppats[k].count391 - gamma_3_391*mpats[k].count391
      end
    end 
  end
  return outp
end



"""
`TermFlag(EmTex::EntireState)`
Takes an entire state (or the empty state for data recording) and returns "true" when every facility has been assigned
the facility for some period.  If everyone has been finished, returns true.
"""
function TermFl(EmTex::EntireState)
  isdone = true
  for mark in keys(EmTex.mkts) # iterates over markets
    isdone = (isdone)&(reduce(&, [ EmTex.mkts[mark].collection[i].finished  for i in keys(EmTex.mkts[mark].collection) ] ))
  end
  return isdone
end


"""
`FindVLBW(demand::Dict, Tex::EntireState, fids::Array{Int64,1})`
This function should take the demand, a state and find every VLBW infant in the market and reassign it
so it will return a Dict with an int for the fid as the key and with the count of VLBW as the value.
We suppose that redistributed VLBW volume is divided equally among the hospitals with high level facilities,
however many there are.  Rounds the volume down to the nearest integer for cases with more than one
facility.
"""
function FindVLBW(demand::Dict{Int64,ProjectModule.LBW}, Tex::EntireState, fids::Array{Int64,1})
  numfac::Int64 = size(fids,1)
  totvlbw::Int64 = 0
  for k1 in keys(Tex.mkts[Tex.fipsdirectory[fids[1]]].collection) # The facilities are always in the same market, so it's ok to take fids[1]
    totvlbw += demand[k1].bt510
    demand[k1].bt510 = 0                                          # Reassign the value here to 0.
    totvlbw += demand[k1].bt1015
    demand[k1].bt1015 = 0
  end
  return outp::Dict{Int64,Float64} = Dict(k => round(Int64,totvlbw/numfac) for k in fids)
end


"""
`MeanCost{T<:Real}(count::T, level::Int; parameters...)`
Should take a patient count T (of transferred patients) and return some quantity of revenue to
be subtracted from the profit of the function.  This will weight these costs according to the
conditional distribution of DRG's in 2005 between 385 and 390 (not 391).  Other assumptions could
be made but this is a start.
"""
function MeanCost{T<:Real}(count::T, level::Int ; alf1::Float64 = 29182.967,
                                                  alf2::Float64 = 22167.6375,
                                                  alf3::Float64 = 23074.8403,
                                                  gamma_1_385::Float64 = 34628.8402,
                                                  gamma_2_385::Float64 = 14921.003,
                                                  gamma_3_385::Float64 = 12822.723,
                                                  gamma_1_386::Float64 = 104578.867,
                                                  gamma_2_386::Float64 = 95366.0004,
                                                  gamma_3_386::Float64 = 69353.471,
                                                  gamma_1_387::Float64 = 34498.5261,
                                                  gamma_2_387::Float64 = 48900.8396,
                                                  gamma_3_387::Float64 = 24639.0552,
                                                  gamma_1_388::Float64 = 26561.8688,
                                                  gamma_2_388::Float64 = 20895.5001,
                                                  gamma_3_388::Float64 = 29775.8381,
                                                  gamma_1_389::Float64 = 20653.5821,
                                                  gamma_2_389::Float64 = 20102.2097,
                                                  gamma_3_389::Float64 = 8279.774,
                                                  gamma_1_390::Float64 = 7372.3301,
                                                  gamma_2_390::Float64 = 2514.8717,
                                                  gamma_3_390::Float64 = 26113.4462,
                                                  weight385::Float64 = 0.065,
                                                  weight386::Float64 = 0.089,
                                                  weight387::Float64 = 0.073,
                                                  weight388::Float64 = 0.161,
                                                  weight389::Float64 = 0.16,
                                                  weight390::Float64 = 0.45)
  if level == 1
    return count*(alf1*count - (weight385*gamma_1_385 + weight386*gamma_1_386 + weight387*gamma_1_387 + weight388*gamma_1_388 + weight389*gamma_1_389 + weight390*gamma_1_390))
  elseif level == 2
     return count*(alf2*count - (weight385*gamma_2_385 + weight386*gamma_2_386 + weight387*gamma_2_387 + weight388*gamma_2_388 + weight389*gamma_2_389 + weight390*gamma_2_390))
  else level == 3
    return count*(alf3*count - (weight385*gamma_3_385 + weight386*gamma_3_386 + weight387*gamma_3_387 + weight388*gamma_3_388 + weight389*gamma_3_389 + weight390*gamma_3_390))
  end
end





"""
`CounterSim(T::Int, Tex::EntireState, pats::patientcollection; lev::Int64 = 1, reassign::Bool = true)`
Runs a counterfactual in the following way: for every hospital in every market, it assigns that hospital a level 3, and assigns every other
hospital a level 1.  Then it runs a simulation of T identical periods - the only variation in the choices is generated by the shock.
The result is a `counterhistory`, which contains for each market a collection of `mkthistory`, each one of which is an array of hospital
fid and a vector of `mktyear` - these contain the distribution of births and the results LBW, VLBW, deaths, etc in the form of a `hyrec`
for each hospital.
By varying the `level` parameter, it is possible to simulate the assignment of level 2 to all hospitals except the special specified one.
Or the assignment of everyone to have level 3.
# TODO: There is at least an interesting counterfactual where everyone has level 3, another where every county does.

Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 20);
patients = NewPatients(Tex);

"""
function CounterSim(T::Int, Tex::EntireState, pats::patientcollection; lev::Int64 = 1, reassign::Bool = true)
  d1 = NewHospDict(Tex) # creates a dict for GenP below.
  d2 = NewHospDict(Tex) # creates a dict for GenM below
  arry1 = zeros(Int64, 1550) # allocates an array for use in GenP.  Can be re-used.
  arry2 = zeros(Int64, 1550) # allocates an array for use in GenM.  Can be re-used.
  res = counterhistory(Dict{Int64, mkthistory}())                                                            # the output - a counterfactual history
  #TODO 02/20/2017 - remove dict comprehension.  
  res.hist = Dict(k => mkthistory(k, Dict{Int64,simrun}()) for k in keys(Tex.mkts))                          # Fill the dictionary with the fids via a comprehension.
  for mk in keys(Tex.mkts) #these are fipscodes
    res.hist[mk].values = Dict( k=>simrun(mk, Dict{Int64,hyrec}(), 0.0, k) for k in keys(Tex.mkts[mk].collection)) # for each FID in the market, a simrun record element.
    for k1 in keys(Tex.mkts[mk].collection)
      res.hist[mk].values[k1].hosprecord = Dict( k2 => hyrec(k2, Array{Int64,1}(T), Array{Int64,1}(T), Array{Int64,1}(T), Array{Float64,1}(T), Array{Float64,1}(T)) for k2 in keys(Tex.mkts[mk].collection))
    end
  end
  CounterCleanResults(res)
  termflag = true                                                                                           # Start the termination flag.
  while termflag
    #TODO - it may make sense to hive this all off and make it a separate function choosing subsets of size n from the hospitals.
    currentfac = Dict{Int64, Int64}()                                                                       # this will be filled with {fipscode, fid} entries for unfinished facilities.
    mkt_fips = Array{Int64,1}()                                                                             # Tracks the fipscodes which still need to be done.
    for el in keys(Tex.mkts)
      if !reduce(&, [Tex.mkts[el].collection[i].finished for i in keys(Tex.mkts[el].collection)])
        pfids = prod(hcat( [ [i, !Tex.mkts[el].collection[i].finished] for i in keys(Tex.mkts[el].collection) ]...) , 1)
        pfid = pfids[findfirst(pfids)]                                                                      # takes the first non-zero element of the above and returns the element.
        currentfac[Tex.fipsdirectory[pfid]] = pfid                                                          # Now the Key is the fipscode and the value is the fid.
        push!(mkt_fips, el)                                                                                 # Want to do this collection of markets which still have uncompleted hospitals.
        for hos in Tex.mkts[el].config
          if hos.fid == pfid
            SetLevel(Tex.mkts[el], pfid, lev)                                                               # Set the level of everyone else in the market to lev, and pfid to 3.
            hos.finished = true
            hos.hasint = true
          else
            hos.hasint = false
          end
        end
      end
    end
    # TODO: stop here - think about cutting off the section above this line below the while termflag component.
    println(mkt_fips)
    UpdateDeterministic(pats)                                                                               # NB: The update happens every time we do a new set of facilities.
    wtpc = WTPMap(pats, Tex)                                                                                # Facilities are unchanging, so WTP will remain constant.
    for i = 1:T                                                                                             # T is now the sim periods, not sequential choices.
      #TODO - Clean these dictionaries up later.    
      GenPChoices(pats, d1, arry1)
      GenMChoices(pats, d2, arry2)
      #TODO - don't reallocate demand every period - clean up 
      mappeddemand, drgp, drgm = PatientDraw(d1, d2, Tex)        # NB: this is creating a Dict{Int64, LBW} of fids and low birth weight volumes.
      pdict = Payoff(drgp, drgm, Tex, wtpc)
      if !reassign                                                                                        # NB: Under this counterfactual, I am restricting investment and NOT transferring.
        for el in mkt_fips
          fac = currentfac[el]
          mortcount = 0.0
          for k in keys(Tex.mkts[el].collection)
            #TODO - this is wrong.  Mappeddemand is NICU admits, NOT all of the demand.
              res.hist[el].values[currentfac[el]].hosprecord[k].totbr[i] = sum(drgp[k])+sum(drgm[k])
              res.hist[el].values[currentfac[el]].hosprecord[k].totlbw[i] = mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510
              res.hist[el].values[currentfac[el]].hosprecord[k].totvlbw[i] = mappeddemand[k].bt1015 + mappeddemand[k].bt510
              res.hist[el].values[currentfac[el]].hosprecord[k].deaths[i] = VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510, Tex.mkts[el].collection[k].level)*(mappeddemand[k].bt1015 + mappeddemand[k].bt510)
              res.hist[el].values[currentfac[el]].hosprecord[k].profit[i] = pdict[k]
              mortcount += res.hist[el].values[currentfac[el]].hosprecord[k].deaths[end]                        # this is confusing - has this behavior changed since 0.4?
          end
          res.hist[el].values[currentfac[el]].yeartot = mortcount
        end
      else     #reassign true                                                                                 #NB: Under this counterfactual, I am restricting investment and *transferring* the VLBW only
        for el in mkt_fips
          fac = currentfac[el]
          mortcount = 0.0
          nofac = Array{Int64,1}()
          hasfac = Array{Int64,1}()
          for fid in keys(Tex.mkts[el].collection)
            if (Tex.mkts[el].collection[fid].level == 1)||(Tex.mkts[el].collection[fid].level == 2) # track which facilities have what?
              push!(nofac, fid)
            else
               push!(hasfac, fid)
            end
          end
          for k in nofac
            res.hist[el].values[currentfac[el]].hosprecord[k].totbr[i] = sum(drgp[k])+sum(drgm[k]) - (mappeddemand[k].bt1015 + mappeddemand[k].bt510)
            res.hist[el].values[currentfac[el]].hosprecord[k].totlbw[i] = mappeddemand[k].bt2025 + mappeddemand[k].bt1520
            res.hist[el].values[currentfac[el]].hosprecord[k].totvlbw[i] = 0
            # TODO: compute the death rate among higher weight babies for fairness, perhaps?
            res.hist[el].values[currentfac[el]].hosprecord[k].deaths[i] = 0.0
            res.hist[el].values[currentfac[el]].hosprecord[k].profit[i] = pdict[k] - MeanCost(mappeddemand[k].bt1015 + mappeddemand[k].bt510,lev)
            mortcount += 0
          end
          for k in hasfac
            sharedvlbw = FindVLBW(mappeddemand, Tex, hasfac)
            res.hist[el].values[currentfac[el]].hosprecord[k].totbr[i] = sum(drgp[k])+sum(drgm[k])
            res.hist[el].values[currentfac[el]].hosprecord[k].totlbw[i] = mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510
            res.hist[el].values[currentfac[el]].hosprecord[k].totvlbw[i] = mappeddemand[k].bt1015 + mappeddemand[k].bt510 + sharedvlbw[k]
            res.hist[el].values[currentfac[el]].hosprecord[k].deaths[i] = VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510 + sharedvlbw[k], Tex.mkts[el].collection[k].level)*(mappeddemand[k].bt1015 + mappeddemand[k].bt510 + sharedvlbw[k])
            res.hist[el].values[currentfac[el]].hosprecord[k].profit[i] = pdict[k] + MeanCost(sharedvlbw[k],3) 
            mortcount += res.hist[el].values[currentfac[el]].hosprecord[k].deaths[end]
          end
          res.hist[el].values[currentfac[el]].yeartot = mortcount
        end
      end
    end
    termflag = !TermFl(Tex)                                                                                      # Will only return "true" when everyone is finished.
  end # of While loop.
  for el in keys(Tex.mkts)
    for k in keys(Tex.mkts[el].collection)
      Tex.mkts[el].collection[k].finished = false                                                                 # Resets all "finished" flags to false.
      Tex.mkts[el].collection[k].level = Tex.mkts[el].collection[k].actuallev                                     # Reassigns to actual level.
    end
  end
  UpdateDeterministic(pats)                                                                                       # Resets the deterministic utility component
  return res
end



"""
`Baseline(T::Int, Tex::EntireState, pats::patientcollection)`
This function runs a baseline simulation of mortality rates in each market.  It simulates for T periods in Order
to reduce the impact that the random component of the choice model has.  It should return a counterhistory for which
a mortality baseline can be determined.
Note: run this on a NEW state record unless the levels have been set back to what they were in reality.
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists);
"""
#TODO: check this against some mortality data from NCHS.
function Baseline(T::Int, Tex::EntireState, pats::patientcollection; levelchange::Bool = false, level::Int64 = 3)
  d1 = NewHospDict(Tex) # creates a dict for GenP below.
  d2 = NewHospDict(Tex) # creates a dict for GenM below
  arry1 = zeros(Int64, 1550) # allocates an array for use in GenP.  Can be re-used.
  arry2 = zeros(Int64, 1550) # allocates an array for use in GenM.  Can be re-used.
  res = counterhistory(Dict{Int64, mkthistory}())                                                        # the output - a counterfactual history
  res.hist = Dict(k => mkthistory(k, Dict{Int64,simrun}()) for k in keys(Tex.mkts))                      # Fill the dictionary with the fids via a comprehension.
  for mk in keys(Tex.mkts)
    res.hist[mk].values = Dict( 0=>simrun(mk, Dict{Int64,hyrec}(), 0.0, 0))                              # 0 records that this is the equilibrium.
    for k1 in keys(Tex.mkts[mk].collection)
      res.hist[mk].values[0].hosprecord = Dict( k2 => hyrec(k2, Array{Int64,1}(T), Array{Int64,1}(T), Array{Int64,1}(T), Array{Float64,1}(T), Array{Float64,1}(T)) for k2 in keys(Tex.mkts[mk].collection))
    end
  end
  CounterCleanResults(res)  
  UpdateDeterministic(pats)                                                                              # NB: The update happens every time we do a new set of facilities.
  wtpc = WTPMap(pats, Tex)                                                                           # Facilities are unchanging, so WTP will remain constant.
  for i = 1:T    
    GenPChoices(pats, d1, arry1)  
    GenMChoices(pats, d2, arry2)                                                                                      # T is now the sim periods, not sequential choices.
    mappeddemand, drgp, drgm = PatientDraw(d1,  d1,  Tex)        # NB: this is creating a Dict{Int64, LBW} of fids and low birth weight volumes.
    pdict = Payoff(drgp, drgm, Tex, wtpc)
    if levelchange                                                                                       # NB: this should permit changing the level of every hospital.
      for k1 in keys(Tex.mkts)
        SetLevel(Tex.mkts[k1], 0, level)                                                                 # NB: this will set the level of *everyone* to "level"
      end
    end
    for el in keys(Tex.mkts)
      fac = 0
      mortcount = 0.0
      for k in keys(Tex.mkts[el].collection)
        res.hist[el].values[0].hosprecord[k].totbr[i] = sum(drgp[k])+sum(drgm[k])
        res.hist[el].values[0].hosprecord[k].totlbw[i] = mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510
        res.hist[el].values[0].hosprecord[k].totvlbw[i] = mappeddemand[k].bt1015 + mappeddemand[k].bt510
        res.hist[el].values[0].hosprecord[k].deaths[i] = VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510, Tex.mkts[el].collection[k].level)*(mappeddemand[k].bt1015 + mappeddemand[k].bt510)
        res.hist[el].values[0].hosprecord[k].profit[i] = pdict[k]
        mortcount += res.hist[el].values[0].hosprecord[k].deaths[end]
      end
      res.hist[el].values[0].yeartot = mortcount
    end
  end
  return res
end






"""
`BaselineCheck(sim::counterhistory)`
This will check the number of VLBW births by market under the baseline simulation.
Care needs to be taken for the counterfactual history due to the different key system in use
in the way the data is recorded.
"""
function BaselineCheck(sim::counterhistory; check = false)
  outp::Dict{Int64, Array{Float64, 1}} = Dict(k => [0.0, 0.0] for k in keys(sim.hist))
  totl::Float64 = 0.0
  var::Float64 = 0.0
  mvars::Array{Float64,1} = Array{Float64,1}()
  for el in keys(sim.hist)
    dvar::Array{Float64,1} = zeros(Float64, size(sim.hist[48453].values[0].hosprecord[4530190].deaths,1))
    for k1 in keys(sim.hist[el].values[0].hosprecord)
      outp[el][1] += mean(sim.hist[el].values[0].hosprecord[k1].totvlbw)
      outp[el][2] += mean(sim.hist[el].values[0].hosprecord[k1].totbr)
      dvar += sim.hist[el].values[0].hosprecord[k1].deaths
      totl += mean(sim.hist[el].values[0].hosprecord[k1].deaths)
    end
    push!(mvars, std(dvar))
  end
  if check
    for k in keys(outp)
      println(k, "  Fraction VLBW ",  outp[k][1]/outp[k][2])
    end
  end
  return outp, totl, mvars
end




"""
`DeathCheck(sim::counterhistory; check = false)`
Per fips code, how many deaths are there among VLBW infants?
This is a check on the validity of the Baseline simulation function.
"""
function DeathCheck(sim::counterhistory; check::Bool = false, houston::Bool = false)
  outp = Dict(k=>0.0 for k in keys(sim.hist))
  for el in keys(sim.hist)
    for k1 in keys(sim.hist[el].values[0].hosprecord)
      outp[el] += mean(sim.hist[el].values[0].hosprecord[k1].deaths)
    end
  end
  if houston #just for debugging print the results of this market.
    vlbw = 0.0
    for k1 in keys(sim.hist[48201].values[0].hosprecord)
      println("---------------------------------------")
      println("Facility ", sim.hist[48201].values[0].hosprecord[k1].fid)
      println("VLBW ", mean(sim.hist[48201].values[0].hosprecord[k1].totvlbw) )
      vlbw += mean(sim.hist[48201].values[0].hosprecord[k1].totvlbw)
      println("Deaths ", mean(sim.hist[48201].values[0].hosprecord[k1].deaths))
      println("Fraction ", mean(sim.hist[48201].values[0].hosprecord[k1].deaths)/mean(sim.hist[48201].values[0].hosprecord[k1].totvlbw))
    end
    println("**********************")
    println("Houston Total: ", outp[48201])
    println("Houston VLBW Total: ", vlbw)
    println("Houston Mortality Rate: ", outp[48201]/vlbw)
  end
  if check
    for k in [48201, 48113, 48453] #keys(outp) # check just the subset.
      println(k, "   Deaths: ", outp[k])
    end
  end
  return outp
end




"""
`BestOutcome(sim::counterhistory, basesim::counterhistory)`
This function looks through all the simulations and finds the one with the lowest
number of deaths relative to the baseline.  The results are a dict with the form
fips => [ Base Deaths, Best Deaths, FID], where FID is the ID of the facility for which
this is best.
"""
function BestOutcome(sim::counterhistory, basesim::counterhistory)
  res = Dict( k=> [0.0, 0.0, 0.0] for k in keys(basesim.hist)) #basedeaths, bestsimdeaths, bestfid.
  for el in keys(basesim.hist)
    res[el][1] = basesim.hist[el].values[0].yeartot
    res[el][2] = basesim.hist[el].values[0].yeartot # Assign this to be the same - see if improvement is possible.
    for k2 in keys(sim.hist[el].values)
      if sim.hist[el].values[k2].yeartot <= res[el][2]
        res[el][2] = sim.hist[el].values[k2].yeartot
        res[el][3] = k2
      end
    end
  end
  return res
end



"""
`SimpleResultsPrint(inp1::Dict{Int64, Array{Float64,1}})`
This doesn't do much - takes the simple version of the results and figures out
how many lives are saved overall in the simple counterfactual.
"""
function SimpleResultsPrint(inp1::Dict{Int64, Array{Float64,1}})
  totsaved::Float64 = 0.0
  noimp = 0
  for k in keys(inp1)
    if inp1[k][3] != 0.0
      totsaved += (inp1[k][1] - inp1[k][2])
    else
      noimp += 1
    end
  end
  return totsaved, noimp
end




"""
`DemandChangeATX(sim::counterhistory, basel::counterhistory, Tex::EntireState)`
Simple function to print some of the outcomes of the equilibrium and non-eq simulations for the
Austin Market - means, stds, etc.
"""
function DemandChangeATX(sim::counterhistory, basel::counterhistory, Tex::EntireState)
  for k in keys(Tex.mkts[48453].collection)
    println("**********************")
    println(Tex.mkts[48453].collection[k].name)
    println("Baseline Mean Admissions: ", mean(basel.hist[48453].values[0].hosprecord[k].totbr))
    println("One NICU Mean Admissions: ", mean(sim.hist[48453].values[k].hosprecord[k].totbr))
    preothermeans = 0.0
    postothermeans = 0.0
    println("Other Firms")
    for k2 in keys(Tex.mkts[48453].collection)
      if k != k2
        println(Tex.mkts[48453].collection[k2].name)
        println("Admissions Pre: ", mean(basel.hist[48453].values[0].hosprecord[k2].totbr))
        println("Admissions Post: ", mean(sim.hist[48453].values[k].hosprecord[k2].totbr))
        preothermeans += mean(basel.hist[48453].values[0].hosprecord[k2].totbr)
        postothermeans += mean(sim.hist[48453].values[k].hosprecord[k2].totbr)
      end
    end
    println("Means of Others PRE: ", preothermeans/6)
    println("Means of Others POST: ", postothermeans/6)
  end
end



"""
`HHI(ch::counterhistory)`
compute the HHI every year.  I should append... the mean and the variance in the HHI for each hospital, probably.
         # lbwcount += ch.hist[k1].values[k2].hosprecord[k3].totlbw[i]
         # vlbwcount += ch.hist[k1].values[k2].hosprecord[k3].totvlbw[i]

"""
function HHI(ch::counterhistory, T::Int64)
  outp::Dict{Int64, Tuple{Float64, Float64}} = Dict{Int64, Tuple{Float64, Float64}}()
  brcount::Int64 = 0
  lbwcount::Int64 = 0
  vlbwcount::Int64 = 0
  hh::Float64 = 0.0
  temparray::Array{Float64,1} = zeros(T) # allocate one array to be rewritten
  for k1 in keys(ch.hist) #iterate over markets.
    for k2 in keys(ch.hist[k1].values) #over hospitals in the market - this is the record when that facility was assigned the level 3
      for i = 1:T # each entry here is a year - I want the year-level.
        hh = 0.0
        brcount = 0
        for k3 in keys(ch.hist[k1].values[k2].hosprecord) # other records in that market in that year when k2 had the level 3
          brcount += ch.hist[k1].values[k2].hosprecord[k3].totbr[i]
        end 
        for k3 in keys(ch.hist[k1].values[k2].hosprecord)
          hh += (ch.hist[k1].values[k2].hosprecord[k3].totbr[i]/brcount)^2
        end 
        temparray[i] = hh 
      end
      outp[k2] = (mean(temparray), std(temparray)) 
      for i = 1:length(temparray)
        temparray[i] = 0.0
      end 
    end 
  end 
  return outp 
end 



"""
`ProfitChange`

- res - counterhistory
- res.hist[FIPS] -> what market are we doing?
- res.hist[FIPS].values[key] → a simrun for the FIPS. Key will be a FID for the hospital with the facility.
- res.hist[FIPS].values[FID] → What FID has the level III ?  
- res.hist[FIPS].values[FID] → This is now a "simrun"
- res.hist[FIPS].values[FID].hosprecord[FID2] → append the results of the simulation for ALL HOSPITALS to the dict of hyrecs in hosprecord
- res.hist[FIPS].values[FID].hosprecord[FID2] → the record for FID2 consists of totbr, totvlbw, deats, profits.


"""

function ProfitChange(ch::counterhistory, T::Int64)
  outp::Dict{Int64, Array{Float64, 1}} = Dict{Int64, Array{Float64, 1}}() # fid, (facility, nofac) ? -> What do I want out of this?
  hasfac::Float64 =0.0
  nofac::Float64 = 0.0
  count::Int64 = 0 # how many times does it appear? TODO - not N-1 for N firms in the market.  Something more than that, right?
  for k1 in keys(ch.hist) # keys of ch are FIPS 
    for k2 in keys(ch.hist[k1].values) # keys of values are FIDs
      for k3 in keys(ch.hist[k1].values[k2].hosprecord) # k3 are fids too. Now I am looking at the hyrecs
        if k3 == ch.hist[k1].values[k2].hasfac  # this is the facility with the fid.
          out[k3][1] += mean(ch.hist[k1].values[k2].hosprecord[k3].profit)
        else 
          out[k3][2] += mean(ch.hist[k1].values[k2].hosprecord[k3].profit) # this one needs to be divided at the end 
        end 
      end 
    end 
  end 
  return outp
end 


"""
`RunCounter1(T::Int64)`
Run the first counterfactual and get the results.
Specifically - runs a Baseline for T periods.  Then compares to the following scenarios:
- Forcing single level 3 with transfers
- Forcing single level 3 w/ out transfers
 BELOW THIS LINE NOT DONE YET: these require adjusted mortality rates.  Fix `VolMortality` in DataStructs.jl
- Assigning everyone to level 3
- Assigning everyone to level 2
- Assigning everyone to level 1
"""
function RunCounter1(T::Int64)
    println("Setting up")
    Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
    CMakeIt(Tex, ProjectModule.fips);
    FillState(Tex, ProjectModule.alldists, 50);
    patients = NewPatients(Tex);
  # Run the baseline then run the counter where each county gets a single level 1 and there *are* transfers.
    println("Running Baseline")
    bl = Baseline(T, Tex, patients)
    println("Running Single Level 3 w/ transfers")
    c1 = CounterSim(T, Tex, patients)
    outc1 = BestOutcome(c1, bl)
    SimpleResultsPrint(outc1)
  # Run the counterfactual where there are no transfers
    Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
    CMakeIt(Tex, ProjectModule.fips);
    FillState(Tex, ProjectModule.alldists, 50);
    patients = NewPatients(Tex);
    println("Running Single level 3 w/ out transfers")
    c1nr = CounterSim(T, Tex, patients; reassign = false)
    outc1nr = BestOutcome(c1nr, bl)
  # # Run the counterfactual where everyone has level 3 (SetLevel fixes the level.)
  #   println("*********")
  #   println("Do not trust these results until the adjusted mortality rate is in use!")
  #   println("Assigning everyone to level 3")
  #   Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
  #   CMakeIt(Tex, ProjectModule.fips);
  #   FillState(Tex, ProjectModule.alldists, 50);
  #   patients = NewPatients(Tex);
  #   c3 = Baseline(T, Tex, patients; levelchange = true, level = 3)
  #   outc3all = BestOutcome(c3, bl)
  # # Run the counterfactual where everyone but the special guy has level 2
  #   println("assigning everyone except the special fac to level 2")
  #   Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
  #   CMakeIt(Tex, ProjectModule.fips);
  #   FillState(Tex, ProjectModule.alldists, 50);
  #   patients = NewPatients(Tex);
  #   c2 = Baseline(T, Tex, patients; levelchange = true, level = 2)
  #   outc2all = BestOutcome(c2, bl)
  # # Run a counterfactual with everyone at level 1
  #   println("Assigning everyone to level 1")
  #   Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
  #   CMakeIt(Tex, ProjectModule.fips);
  #   FillState(Tex, ProjectModule.alldists, 50);
  #   patients = NewPatients(Tex);
  #   c1all = Baseline(T, Tex, patients; levelchange = true, level = 1)
  #   outc1all = BestOutcome(c1all, bl)
  return bl, c1, c1nr, outc1, outc1nr, outc3all, outc2all, outc1all
end



 baseline, c1, c1nr, outc1, outc1nr, outc3all, outc2all, outc1all = RunCounter1(5)


#=

Documentation of Record Types for Counter 1:

- Top level: RunCounter1 returns three "counterhistory" objects.  Call this "CH".
- Counterhistory CH has one field: a dict called "hist"
- CH.hist[key] to access. 
- The keys of "CH.hist" are FIPS codes.
- Each CH.hist[key] is a "mkthistory" type object.
- mkthistory MH has fields fips::Int64 and values::Dict{Int64, simrun}
- MH.fips::Int64 is a fipscode
- MH.values[key] is a "simrun", where MH is CH.hist[fipscode]
- CH.hist[FIPS key].values[FID key] returns the "simrun" SR
- Each "simrun" SR is a collection: a fips code, hosprecord dict{Int64, hyrec}, yeartot, and hasfac (an Int64)
- SR.fips::Int64 is a fipscode, 
- SR.hosprecord::Dict{Int64, hyrec} is a collection of hyrecs (hospital-year records)
- SR.yeartot::Float64 records total mortality in the year 
- SR.hasfac::Int64 is a fid
- hyrec HY is a collection of: a fid, totbr total births, totlbw total low birth weight, totvlbw total low birth weight, deaths and profit
- HY.fid::Array{Int64,1}
- HY.totbr::Array{Int64,1}
- HY.totvlbw::Array{Int64,1}
- HY.deaths::Array{Int64,1}
- HY.profit::Array{Float64,1}
- Or this can be accessed with HY replaced with: CH.hist[FIPSCODE].values[FID].hosprecord[FID]
- What happens in a simulation?  See the CounterSim function in the section starting "from i = 1:T"
- Two cases: transferring and no transferring.  Under transferring all VLBW are counted and sent to hospital with facility
- Else every hospital gets a record.  How stored?  
- res.hist[FID].values[currentfac[el]].hosprecord[k] 
- currentfac[el] is the FID of a hospital w/ the facility.  
- k is the fid of the hospital whose record is being created (which may or may not be the hospital with the facility)
- res - counterhistory
- res.hist[FIPS] -> what market are we doing?
- res.hist[FIPS].values[key] → a simrun for the FIPS. Key will be a FID for the hospital with the facility.
- res.hist[FIPS].values[FID] → What FID has the level III ?  
- res.hist[FIPS].values[FID] → This is now a "simrun"
- res.hist[FIPS].values[FID].hosprecord[FID2] → append the results of the simulation for ALL HOSPITALS to the dict of hyrecs in hosprecord
- res.hist[FIPS].values[FID].hosprecord[FID2] → the record for FID2 consists of totbr, totvlbw, deats, profits.

TODO - getting Δ profits, costs, HHI depends on these comparisons.  But now I know how they can be done.  








Do what with that?
- Run BaselineCheck on baseline: (dictionary of deaths per county , total deaths, variances across markets) = BaselineCheck(baseline)
-

=#






#=


### NB: UNCOMPLETED .... Below this line.

"""
`MortalityGet(baseline::counterhistory)`
- Computes the mortality within markets and in the whole state for each of the two sims.
- What do I want this to return?
- For each market, a list of the unique keys.  That's a simulation run.
- Also for each market, how many hospitals are there?
- Return the average number of deaths in the market.
- Really I only care about the best one.
"""
#TODO - this is going to require some kind of record type, I'll bet.  The potential dimension of the return is really variable
# especially if I start selecting 2 or 3 or however many hospitals to have level 3.
# The way this gets filled up is going to change depending on exactly what counterfactual is being run here.
function MortalityGet(sim::counterhistory, base::counterhistory)
  outp = StateResults(Dict{Int64,MktResultsReport}(), Dict{Int64, MktResultsReport}())
  outp.mkts = Dict(k=>MktResultsReport(k, Array{Int64,1}(), 0.0, Dict{Int64, HResults}() ) for k in keys(sim.hist))
  for k in keys(outp.mkts)
    outp.mkts[k].simres = Dict(k1=> HResults(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,k1,0,false) for k1 in keys(sim.hist[k].values))
  end
  outp.baseln = Dict(k=>MktResultsReport(k, Array{Int64, 1}(), 0.0, Dict{Int64, HResults}() ) for k in keys(base.hist))
  for k in keys(base.history)
    outp.baseln[k].simres = Dict(k1=>HResults(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,k1,0,false) for k1 in keys(base.hist[k].values))
  end
  for

  end
  return outp
end



type HResults
  brmean::Float64
  brstd::Float64
  lbwm::Float64
  lbwstd::Float64
  vlbwm::Float64
  vlbwstd::Float64
  prm::Float64
  prstd::Float64
  fid::Int64
  level::Int64
  baseline::Bool
end




type MktResultsReport
  fips::Int64
  hadinten::Array{Float64,1}
  deaths::Float64
  simres::Dict{Int64, HResults}
end

type StateResults
  mkts::Dict{Int64, MktResultsReport}
  baseln::Dict{Int64, MktResultsReport}
end



### To be run:
#=
1.  Everyone at level 1.
2.  Everyone at level 2.
3.  Everyone at level 3.
4.  N-1 at level 1, one at level 3.  (For markets with more than 2 hospitals.)
  - With transfers
  - Without transfers
5.  N-2 at level 1, one at level 2, one at level 3.  No transfers. (For markets with more than 3 hospitals.)
6.  N-1 at level 2, one at level 3.
  - With transfers
  - Without transfers
7.  The baseline, actual scenario.
8.  The best two hospitals at level 3.  No transfers.
9.  Regionalized - Set up regionalized markets, transfers allowed - ???
=#







# These may be too big - or they may be fine.  Check them against the means and variances of the medicaid reimbursements for the same years.
α₁  29182.967
α₂  22167.6375
α₃  23074.8403

# 385: Neonates - died or transferred to another facility.  -> 789: Neonates died or transferred.
γ¹₅  34628.8402
γ₅²  14921.003
γ₅³  12822.723

# 386: Extreme immaturity or resp. distress syndrome -> 790: Extreme immaturity or resp. distress syndrome.
γ₆¹  104578.867
γ₆²  95366.0004
γ₆³  69353.471

# 387: Prematurity w/ major problems -> 791: Prematurity w/ major problems.
γ₇¹  34498.5261
γ₇²  48900.8396
γ₇³  24639.0552

# 388: Prematurity w/out major problems -> 792: Prematurity w/out major problems.
γ₈¹  26561.8688
γ₈²  20895.5001
γ₈³  29775.8381

# 389: Full Term Neonate w/ major problems -> 793: Full Term Neonate w/ significant Problems
γ₉¹  20653.5821
γ₉²  20102.2097
γ₉³  8279.774

# 390: Neonate Significant Problems -> 794: Neonate w/ other significant problems.
γ₀¹  7372.3301
γ₀²  2514.8717
γ₀³  26113.4462

# 391: Normal Newborn.  -> 795: Normal neonate, bwt > 2499 g.
γ₁¹  27018.9915
γ₁²  15079.2889
γ₁³  1912.7285

ϕ12  115855.8434
ϕ13  216801.3608
ϕ1EX  -6509.848
ϕ21  123793.385
ϕ23  235831.3169
ϕ2EX  -14860.5887
ϕ31  85109.5513
ϕ32  189956.7601
ϕ3EX  -38939.4797
=#
