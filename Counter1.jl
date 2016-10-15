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
Computes the actual firm payoffs.
"""
function Payoff(hos::hospital; params = [] ) # params


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


function CounterSim(T::Int, Tex::EntireState, pats::patientcollection)
  # Runs a T period simulation using the whole state and whole collection of patient records.
  # It's easy enough to run the sim 20 times for each hospital as the one with the NICU.
  # TODO: There is at least an interesting counterfactual where everyone has level 3, another where every county does.
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
            SetLevel(Tex.mkts[el], pfid)                                                         # Set the level of everyone else in the market to 1, and pfid to 3.
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
    for i = 1:T                                                                                 # T is now the sim periods, not sequential choices.
      mappeddemand = PatientDraw(GenPChoices(pats, Tex),  GenMChoices(pats, Tex),  Tex )        # NB: this is creating a Dict{Int64, LBW} of fids and low birth weight volumes.
      for el in mkt_fips
        fac = currentfac[el]
        myr = mktyear(el, Dict{Int64, hyrec}(), 0, fac)
        for k in keys(Tex.mkts[el].collection)
          myr.hosprecord[k] = hyrec(k,
                                    sum(mappeddemand[k]),
                                    mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                    mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                    floor(VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510, Tex.mkts[el].collection[k].level)*(mappeddemand[k].bt1015 + mappeddemand[k].bt510)))
          if Tex.mkts[el].collection[k].hasint  # test that hospital is the one w/ fac.
            myr.hasfac = k
          end
        end
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
  res = counterhistory(Dict{Int64, mkthistory}())                                               # the output - a counterfactual history
  res.hist = Dict(k => mkthistory(k, Dict{Int64,Array{mktyear,1}}()) for k in keys(Tex.mkts))   # Fill the dictionary with the fids via a comprehension.
# #TODO: get rid of this - there are not N sims for the N hospitals anymore.
  for mk in keys(Tex.mkts)
    res.hist[mk].history = Dict( 0=>Array{mktyear,1}())                                       # 0 records that this is the equilibrium.
  end
  UpdateDeterministic(pats)                                                                   # NB: The update happens every time we do a new set of facilities.
  for i = 1:T                                                                                 # T is now the sim periods, not sequential choices.
    mappeddemand = PatientDraw(GenPChoices(pats, Tex),  GenMChoices(pats, Tex),  Tex )        # NB: this is creating a Dict{Int64, LBW} of fids and low birth weight volumes.
    for el in keys(Tex.mkts)
      fac = 0 # this has to be changed!
      myr = mktyear(el, Dict{Int64, hyrec}(), 0, fac)
      for k in keys(Tex.mkts[el].collection)
        myr.hosprecord[k] = hyrec(k,
                                  sum(mappeddemand[k]),
                                  mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                  mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                  floor(VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510, Tex.mkts[el].collection[k].level)*(mappeddemand[k].bt1015 + mappeddemand[k].bt510)))
      end
      push!(res.hist[el].history[fac], myr)                                                              # NB: at this point, we have a market-year record with each hospital recorded.  Add it to the market history within the counterhistory
    end
  end
  return res
end




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
