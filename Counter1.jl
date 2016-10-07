# First Counterfactual.

#=

Outline - track demand,
track fraction LBW,
VLBW,
Actually don't really care about demand... only VLBW, LBW
What DRG's?
All of 385, 386, 387
what about 388?
All of 389  390
None of 391
=#

include("DataStructs.jl")


# Define an object which is a fid, total volume and total mortality.
# why shouldn't I just add these to datastructs?
# I don't need them.  Can just define new ones.
# well, think about what I need first.
# I can define a larger abstract type and then hospital and chospital as subtypes?



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
`FillState(Tex::EntireState)`
fills the entire state record with elements of the chospital type
for the counterfactual only.
Note that this needs to be called AFTER the function `CMakeIt(Tex::EntireState, fip::Vector)` is called
on an empty state record.
This version also adds all of the `chospitals` directly to the `Market.collection` dictionary.
"""
function FillState(Tex::EntireState, data::Matrix; lev105loc = 97, lev205loc = 98, lev305loc = 99, lev1515loc = 101, lev2515loc = 102, lev3515loc = 103, lev11525loc = 105, lev21525loc = 106, lev31525loc = 107)
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
      ProjectModule.chospital( data[i, 74],
                data[i,94],
                data[i, 95],
                data[i, 82],
                fips,
                level,
                Array{Int64,1}(), #volume
                Array{Int64, 1}(), #mortality
                Array{Float64,1}(), #ppayoff
                Array{Float64,1}(), #mpayoff
                  0    , # beds added later.
                LBW(Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}(), Array{Int64,1}()) # LBW Infants.  
                false ) )
      Tex.mkts[fips].collection[data[i,74]] = ProjectModule.chospital( data[i, 74],
                data[i,94],
                data[i, 95],
                data[i, 82],
                fips,
                level,
                Array{Int64,1}(), #volume
                Array{Int64, 1}(), #mortality
                Array{Float64,1}(), #ppayoff
                Array{Float64,1}(), #mpayoff
                  0    , # beds added later.
                false )
    end
    # push all hospital fid/ fips pairs into the directory.
    Tex.fipsdirectory[data[i, 74]] = fips # now for the whole state I can immediately figure out which market a hospital is in.
  end
  return Tex
end

# Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
# CMakeIt(Tex, ProjectModule.fips);
# FillState(Tex, ProjectModule.data05)


"""
`Mortality(mkt::Market)`
Computes the mortality rate at the market level.
"""
function Mortality(mkt::Market; prob)


end

"""
`Payoff(hos::hospital; params = [])`
Computes the actual firm payoffs.
"""
function Payoff(hos::hospital; params = [] ) # params


end


function CounterSim(T::Int, Tex::EntireState, pats::patientcollection; entrants = [0, 1], entryprobs = [0.9895, 0.0105] )
  # Runs a T period simulation using the whole state and whole collection of patient records.
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

    end
    UpdateDeterministic(pats)                                                                  # Updates deterministic component of utility
  end
  return Tex                                                                                   # Returns the whole state so the results can be written out.
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
