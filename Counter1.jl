# First Counterfactual.





# Create and fill state:
# Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
# CMakeIt(Tex, ProjectModule.fips);
# FillState(Tex, ProjectModule.data05);


# Patient collection:
# patients = NewPatients(Tex);



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


function CounterSim(T::Int, Tex::EntireState, pats::patientcollection; entrants = [0, 1], entryprobs = [0.9895, 0.0105] )
  # Runs a T period simulation using the whole state and whole collection of patient records.
  # It's easy enough to run the sim 20 times for each hospital as the one with the NICU.
  res = counterhistory(Dict{Int64, mkthistory}())                                               # the output - a counterfactual history
  # NB: switch the order here - do 20 years of each market, but maybe a subset of markets only.
  for el in Tex.ms
    #TODO: make some arrangement to make sure there is one equilibrium simulation here
    if size(el.config,1)>1 # NB: There is at least an interesting counterfactual where everyone has level 3, another where every county does.
      mkh = mkthistory(el.fipscode, Dict{Int64, mktyear}())  #NB: The dict is indexed by fids.
      unf = FindUndone(el) # NB: collection of unfinished fids.
      if size(unf,1)>1
        ufid = unf[1] # Pick the first undone fid
        SetLevel(el, ufid)
        for i = 1:T # T is now the sim periods, not sequential choices.
          myr = mktyear(el.fipscode, Dict{Int64, hyrec}(), 0)           # Create an empty market-year record. Dict contains fid/patient volumes.
          UpdateDeterministic(pats)
          #TODO: PDemandMap and MDemandMap are trying to map to fields which chospital doesn't have.
          mappeddemand = PatientDraw(GenPChoices(pats, Tex),  GenMChoices(pats, Tex),  Tex ) #NB: this is creating a Dict{Int64, LBW} of fids and low birth weight volumes.
          for k in keys(mappeddemand)
            myr.hosprecord[k] = hyrec(k,
                                      sum(mappeddemand[k]),
                                      mappeddemand[k].bt2025 + mappeddemand[k].bt1520 + mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                      mappeddemand[k].bt1015 + mappeddemand[k].bt510,
                                      floor(VolMortality(mappeddemand[k].bt1015 + mappeddemand[k].bt510, el.collection[k])*(mappeddemand[k].bt1015 + mappeddemand[k].bt510)))
          end
          entrant = sample(entrants, WeightVec(entryprobs))
          #TODO: the entrants have to be removed this time.
          if entrant != 0
            entloc = NewEntrantLocation(el)                                                        # called on the market
            newfid = -floor(rand()*1e6)-1000000                                                    # all entrant fids negative to facilitate their removal later.
            entr = ProjectModule.chospital( data[i, 74],
                      data[i,94],
                      data[i, 95],
                      data[i, 82],
                      fips,
                      level,
                      Array{Int64,1}(), #volume
                      Array{Int64, 1}(), #mortality
                      Array{Float64,1}(), #ppayoff
                      Array{Float64,1}(), #mpayoff
                        0    , # TODO: beds added later.
                      LBW(0,0,0,0,0,0), # LBW Infants.
                      false, # has intensive
                      false )
            push!(el.config, entr)                                                                 # need to create a new record for this hospital in the market
            el.collection[newfid] = entr
            for elm in el.config                                                                   # need to add it to the dictionary too:
              NeighborAppend(elm, entr)
              NeighborAppend(entr, elm)
            end
             HospUpdate(entr, entrant)
             mkh.history[i] = mktyear                                                          #"entrant" is the level
          end
        end
      end
    end                                                                 # Updates deterministic component of utility
  end
  return res                                                                                  # Returns the collection of year results.
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
