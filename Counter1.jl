# First Counterfactual.





# Create and fill state and Patient collection:
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.data05);
patients = NewPatients(Tex);






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
`Payoff(hos::hospital; params = [])`
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
  outp = Dict(k => 0.0 for k in keys(ppats))
  for k in keys(ppats) # multiplied by patient volumes or not?
    if Tex.mkts[Tex.fipsdirectory[k]].collection[k].level == 1
      outp[k] = alf1*wtp[k]*(sum(ppats[k])+sum(mpats[k])) - gamma_1_385*ppats[k].count385 - gamma_1_385*mpats[k].count385 - gamma_1_386*ppats[k].count386 - gamma_1_386*mpats[k].count386 - gamma_1_387*ppats[k].count387 - gamma_1_387*mpats[k].count387 - gamma_1_388*mpats[k].count388 - gamma_1_388*ppats[k].count388 - gamma_1_389*mpats[k].count389 - gamma_1_389*ppats[k].count389 - gamma_1_390*ppats[k].count390 - gamma_1_390*mpats[k].count390 - gamma_1_391*ppats[k].count391 - gamma_1_391*mpats[k].count391
    elseif Tex.mkts[Tex.fipsdirectory[k]].collection[k].level == 2
      outp[k] = alf2*wtp[k]*(sum(ppats[k])+sum(mpats[k])) - gamma_2_385*ppats[k].count385 - gamma_2_385*mpats[k].count385 - gamma_2_386*ppats[k].count386 - gamma_2_386*mpats[k].count386 - gamma_2_387*ppats[k].count387 - gamma_2_387*mpats[k].count387 - gamma_2_388*mpats[k].count388 - gamma_2_388*ppats[k].count388 - gamma_2_389*mpats[k].count389 - gamma_2_389*ppats[k].count389 - gamma_2_390*ppats[k].count390 - gamma_2_390*mpats[k].count390 - gamma_2_391*ppats[k].count391 - gamma_2_391*mpats[k].count391
    else # level is 3
      outp[k] = alf3*wtp[k]*(sum(ppats[k])+sum(mpats[k])) - gamma_3_385*ppats[k].count385 - gamma_3_385*mpats[k].count385 - gamma_3_386*ppats[k].count386 - gamma_3_386*mpats[k].count386 - gamma_3_387*ppats[k].count387 - gamma_3_387*mpats[k].count387 - gamma_3_388*mpats[k].count388 - gamma_3_388*ppats[k].count388 - gamma_3_389*mpats[k].count389 - gamma_3_389*ppats[k].count389 - gamma_3_390*ppats[k].count390 - gamma_3_390*mpats[k].count390 - gamma_3_391*ppats[k].count391 - gamma_3_391*mpats[k].count391
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
`CounterSim(T::Int, Tex::EntireState, pats::patientcollection; lev::Int64 = 1)`
Runs a counterfactual in the following way: for every hospital in every market, it assigns that hospital a level 3, and assigns every other
hospital a level 1.  Then it runs a simulation of T identical periods - the only variation in the choices is generated by the shock.
The result is a `counterhistory`, which contains for each market a collection of `mkthistory`, each one of which is an array of hospital
fid and a vector of `mktyear` - these contain the distribution of births and the results LBW, VLBW, deaths, etc in the form of a `hyrec`
for each hospital.
"""
function CounterSim(T::Int, Tex::EntireState, pats::patientcollection; lev::Int64 = 1)
  # Runs a T period simulation using the whole state and whole collection of patient records.
  # It's easy enough to run the sim 20 times for each hospital as the one with the NICU.
  # TODO: There is at least an interesting counterfactual where everyone has level 3, another where every county does.
  # TODO: and the main counterfactual can be varied where the one special hosp has level 3 and everyone else has level 1, vs. the same where everyone else has level 2.
  res = counterhistory(Dict{Int64, mkthistory}())                                               # the output - a counterfactual history
  res.hist = Dict(k => mkthistory(k, Dict{Int64,Array{mktyear,1}}()) for k in keys(Tex.mkts))   # Fill the dictionary with the fids via a comprehension.
  for mk in keys(Tex.mkts)
    res.hist[mk].history = Dict( k=>Array{mktyear,1}() for k in keys(Tex.mkts[mk].collection))
  end
  termflag = true                                                                               # Start the termination flag.
  while termflag
    currentfac = Dict{Int64, Int64}()                                                           # this will be filled with {fipscode, fid} entries for unfinished facilities.
    mkt_fips = Array{Int64,1}()                                                                 # Tracks the fipscodes which still need to be done.
    for el in keys(Tex.mkts)
      if !reduce(&, [Tex.mkts[el].collection[i].finished for i in keys(Tex.mkts[el].collection)])
        pfids = prod(hcat( [ [i, !Tex.mkts[el].collection[i].finished] for i in keys(Tex.mkts[el].collection) ]...) , 1)
        pfid = pfids[findfirst(pfids)]                                                           # takes the first non-zero element of the above and returns the element.
        currentfac[Tex.fipsdirectory[pfid]] = pfid                                               # Now the Key is the fipscode and the value is the fid.
        push!(mkt_fips, el)                                                                  # Want to do this collection of markets which still have uncompleted hospitals.
        for hos in Tex.mkts[el].config
          if hos.fid == pfid
            SetLevel(Tex.mkts[el], pfid, lev)                                                         # Set the level of everyone else in the market to lev, and pfid to 3.
            hos.finished = true
            hos.hasint = true
          else
            hos.hasint = false
          end
        end
      end
    end
    println(mkt_fips)
    UpdateDeterministic(pats)                                                                   # NB: The update happens every time we do a new set of facilities.
    wtpc = WTPMap(patients, Tex)                                                                # Facilities are unchanging, so WTP will remain constant.
    for i = 1:T                                                                                 # T is now the sim periods, not sequential choices.
      mappeddemand, drgp, drgm = PatientDraw(GenPChoices(pats, Tex),  GenMChoices(pats, Tex),  Tex )        # NB: this is creating a Dict{Int64, LBW} of fids and low birth weight volumes.
      pdict = Payoff(drgp, drgm, Tex, wtpc)
      for el in mkt_fips
        fac = currentfac[el]
        myr = mktyear(el, Dict{Int64, hyrec}(), 0.0, fac)
        mortcount = 0.0
        for k in keys(Tex.mkts[el].collection)
          myr.hosprecord[k] = hyrec(k,
                                    sum(mappeddemand[k]),
                                    mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                    mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                    VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510, Tex.mkts[el].collection[k].level)*(mappeddemand[k].bt1015 + mappeddemand[k].bt510),
                                    pdict[k])
          mortcount += myr.hosprecord[k].deaths
          if Tex.mkts[el].collection[k].hasint  # test that hospital is the one w/ fac.
            myr.hasfac = k
          end
        end
        myr.yeartot = mortcount
        push!(res.hist[el].history[fac], myr)                                                              # NB: at this point, we have a market-year record with each hospital recorded.  Add it to the market history within the counterhistory
      end
    end
    termflag = !TermFl(Tex)                                                                    # Will only return "true" when everyone is finished.
  end
  return res
end



"""
`Baseline(T::Int, Tex::EntireState, pats::patientcollection)`
This function runs a baseline simulation of mortality rates in each market.  It simulates for T periods in Order
to reduce the impact that the random component of the choice model has.  It should return a counterhistory for which
a mortality baseline can be determined.
"""
#TODO: check this against some mortality data from NCHS.
function Baseline(T::Int, Tex::EntireState, pats::patientcollection)
  res = counterhistory(Dict{Int64, mkthistory}())                                                        # the output - a counterfactual history
  res.hist = Dict(k => mkthistory(k, Dict{Int64,Array{mktyear,1}}()) for k in keys(Tex.mkts))            # Fill the dictionary with the fids via a comprehension.
  for mk in keys(Tex.mkts)
    res.hist[mk].history = Dict( 0=>Array{mktyear,1}())                                                  # 0 records that this is the equilibrium.
  end
  UpdateDeterministic(pats)                                                                              # NB: The update happens every time we do a new set of facilities.
  wtpc = WTPMap(patients, Tex)                                                                           # Facilities are unchanging, so WTP will remain constant.
  for i = 1:T                                                                                            # T is now the sim periods, not sequential choices.
    mappeddemand, drgp, drgm = PatientDraw(GenPChoices(pats, Tex),  GenMChoices(pats, Tex),  Tex)        # NB: this is creating a Dict{Int64, LBW} of fids and low birth weight volumes.
    pdict = Payoff(drgp, drgm, Tex, wtpc)
    for el in keys(Tex.mkts)
      fac = 0 # this has to be changed!
      myr = mktyear(el, Dict{Int64, hyrec}(), 0.0, fac)
      mortcount = 0.0
      for k in keys(Tex.mkts[el].collection)
        myr.hosprecord[k] = hyrec(k,
                                  sum(mappeddemand[k]),
                                  mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                  mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                  VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510, Tex.mkts[el].collection[k].level)*(mappeddemand[k].bt1015 + mappeddemand[k].bt510),
                                  pdict[k])
        mortcount += myr.hosprecord[k].deaths
      end
      myr.yeartot = mortcount
      push!(res.hist[el].history[fac], myr)                                                              # NB: at this point, we have a market-year record with each hospital recorded.  Add it to the market history within the counterhistory
    end
  end
  return res
end



"""
`MortalityGet(baseline::counterhistory)`
- Computes the mortality within markets and in the whole state for each of the two sims.
- What do I want this to return?
- For each market, a list of the unique keys.  That's a simulation run.
- Also for each market, how many hospitals are there?
- Return the average number of deaths in the market.

"""
#TODO - this is going to require some kind of record type, I'll bet.  The potential dimension of the return is really variable
# especially if I start selecting 2 or 3 or however many hospitals to have level 3.
function MortalityGet(sim::counterhistory; simlen = 20)
  outp =
  for k1 in keys(sim.hist)
    for k2 in keys(sim.hist[k1].history)

    end
  end

end

type HResults
  bmean::Float64
  bstd::Float64
  lbwm::Float64
  lbwstd::Float64
  vlbwm::Float64
  vlbwstd::Float64
  pm::Float64
  pstd::Float64
  fid::Int64
  level::Int64
  baseline::Bool 
end


type MktResultsReport
  fips::Int64
  hadinten::Array{Float64,1}
  deaths::Float64
  res::Dict{Int64, HResults}
end





# NB: these numbers give magnitudes which are way too big.

α₁  29182.967
α₂  22167.6375
α₃  23074.8403
γ¹₅  34628.8402
γ₅²  14921.003
γ₅³  12822.723
γ₆¹  104578.867
γ₆²  95366.0004
γ₆³  69353.471
γ₇¹  34498.5261
γ₇²  48900.8396
γ₇³  24639.0552
γ₈¹  26561.8688
γ₈²  20895.5001
γ₈³  29775.8381
γ₉¹  20653.5821
γ₉²  20102.2097
γ₉³  8279.774
γ₀¹  7372.3301
γ₀²  2514.8717
γ₀³  26113.4462
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
