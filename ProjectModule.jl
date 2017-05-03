
#lis = addprocs(1)

module ProjectModule
  #using DataFrames # this needs to be inside and outside the module?
  using Distributions
  using StatsBase
  using Combinatorics

  dir = pwd()
  global pathdata = "";global pathpeople = "";global  pathprograms = "";
  if (dir == "/Users/austinbean/Desktop/dynhosp")||(dir=="/Users/austinbean/julia6dev/julia")
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
    global pathdata = "/home1/04179/abean/dynhosp/"
    global pathpeople = "/home1/04179/abean/dynhosp/"
    global pathprograms = "/home1/04179/abean/dynhosp/"
  elseif (dir=="/work/04179/abean/dynhosp")
    global pathdata = "/work/04179/abean/dynhosp/"
    global pathpeople = "/work/04179/abean/dynhosp/"
    global pathprogram = "/work/04179/abean/dynhosp/"
  else
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
    println("Hey you're in the wrong directory!")
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  end

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

  type LBW
    bt05::Int64
    bt510::Int64
    bt1015::Int64
    bt1520::Int64
    bt2025::Int64
    bt2580::Int64
  end

  import Base.sum
  function sum(x::LBW)
    return x.bt05 + x.bt510 + x.bt1015 + x.bt1520 + x.bt2025 + x.bt2580
  end


  type neighbors
    level105::Int64
    level205::Int64
    level305::Int64
    level1515::Int64
    level2515::Int64
    level3515::Int64
    level11525::Int64
    level21525::Int64
    level31525::Int64
  end

 Base.start(::neighbors) = :level105

function Base.next(n::neighbors, state)
  if state == :level105
    return n.level105, :level205
  elseif state == :level205
    return n.level205, :level305
  elseif state == :level305
    return n.level305, :level1515
  elseif state == :level1515
    return n.level1515, :level2515
  elseif state == :level2515
    return n.level2515, :level3515
  elseif state == :level3515
    return n.level3515, :level15125
  elseif state == :level15125
    return n.level11525, :level21525
  elseif state == :level21525
    return n.level21525, :level31525
  elseif state == :level31525
    return n.level31525, :level4
  else
    return :level4
  end
end

function Base.done(n::neighbors, state)
  if state == :level4
    return true
  else
    return false
  end
end







  function Base.isequal(n1::neighbors, n2::neighbors)::Bool
    n1.level105 == n2.level105 && n1.level205 == n2.level205 && n1.level305 == n2.level305 && n1.level1515 == n2.level1515 && n1.level2515 == n2.level2515 && n1.level3515 == n2.level3515 && n1.level11525 == n2.level11525 && n1.level21525 == n2.level21525 && n1.level31525 == n2.level31525
  end

  function Base.hash(n1::neighbors)
    hash((n1.level105, n1.level205, n1.level305, n1.level1515, n1.level2515, n1.level3515, n1.level11525, n1.level21525, n1.level31525) )
  end

  function Base.sum(n::neighbors)
    n.level105+n.level205+n.level305+n.level1515+n.level2515+n.level3515+n.level11525+n.level21525+n.level31525
  end
#TODO - extend to take arbitrary args.  This should require splatting I think.
  function Base.sum(n1::neighbors, n2::neighbors)
    neighbors(n1.level105+n2.level105, n1.level205+n2.level205, n1.level305+n2.level305, n1.level1515+n2.level1515, n1.level2515+n2.level2515, n1.level3515+n2.level3515, n1.level11525+n2.level11525, n1.level21525+n2.level21525, n1.level31525+n2.level31525)
  end

  abstract type Fac end # 6devfix

  immutable initial
    level::Int64
  end


  type hospital <: Fac
    fid::Int64
    lat::Float64
    long::Float64
    name::AbstractString
    fipscode::Int64
    level::Int64
    init::initial # initial level is immutable.
    levelhistory::Array{Int64,1}
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

# this hospital type is for the counterfactual only.
#TODO - add a field for the actual level, so it can be reset.
  type chospital <: Fac
    fid::Int64
    lat::Float64
    long::Float64
    name::AbstractString
    fipscode::Int64
    level::Int64
    actuallev::Int64
    totalv::Array{Int64}
    mortality::Array{Int64,1}
    ppayoff::Array{Float64,1}
    mpayoff::Array{Float64,1}
    bedcount::Float64
    lbinf::LBW
    hasint::Bool
    finished::Bool
  end


  type Market{T<:Fac}
    config::Array{T, 1}
    collection::Dict{Int64, T} # create the dict with a comprehension to initialize
    fipscode::Int64
    noneqrecord::Dict{Int64, Bool}
  end


  type EntireState
    ms::Array{Market, 1}
    mkts::Dict{Int64, Market}   # Link markets by FIPS code via dictionary.
    fipsdirectory::Dict{Int64,Int64} # Directory should be hospital fid / market fips
  end

          #### NB: Demand-side Data Structures ######

  mutable struct patientcount{T<:Real}
   count385::T
   count386::T
   count387::T
   count388::T
   count389::T
   count390::T
   count391::T
  end

  import Base.+
  function +(x::patientcount, y::patientcount)
    return patientcount(x.count385 + y.count385, x.count386 + y.count386, x.count387 + y.count387, x.count388 + y.count388, x.count389 + y.count389, x.count390 + y.count390, x.count391 + y.count391)
  end

  function sum(x::patientcount)
    return x.count385 + x.count386 + x.count387 + x.count388 + x.count389 + x.count390 + x.count391
  end

  Base.start(::patientcount) = :count385
  
  function Base.next(p::patientcount, state)
    if state == :count385
      return p.count385, :count386
    elseif state == :count386
      return p.count386, :count387
    elseif state == :count387
      return p.count387, :count388
    elseif state == :count388
      return p.count388, :count389
    elseif state == :count389
      return p.count389, :count390
    elseif state == :count390
      return p.count390, :count391
    elseif state == :count391
      return p.count391, :count392
    else
      return :count392
    end
  end
  function Base.done(p::patientcount, state)
    if state == :count392
      return true
    else
      return false
    end
  end

  mutable struct patientrange{T<:Real}
    l385::T 
    u385::T
    l386::T 
    u386::T    
    l387::T 
    u387::T    
    l388::T 
    u388::T    
    l389::T 
    u389::T    
    l390::T 
    u390::T    
    l391::T 
    u391::T
  end 

  import Base.+
  function +(p1::patientrange, p2::patientrange)
    return patientrange(p1.l385 + p2.l385, p1.u385 + p2.u385, p1.l386 + p2.l386, p1.u386 + p2.l386, p1.l387 + p2.l387, p1.u387 + p2.l387, p1.l388 + p2.l388, p1.u388 + p2.l388, p1.l389 + p2.l389, p1.u389 + p2.l389, p1.l390 + p2.l390, p1.u390 + p2.l390, p1.l391 + p2.l391, p1.u391 + p2.l391)
  end 

  function DrawPatients(p::patientrange, i::Int64)::Int64
    if i == 385
      return rand(p.l385:p.u385)
    elseif i == 386
      return rand(p.l386:p.u386)
    elseif i == 387
      return rand(p.l387:p.u387)
    elseif i == 388
      return rand(p.l388:p.u388)
    elseif i == 389
      return rand(p.l389:p.u389)
    elseif i == 390
      return rand(p.l390:p.u390)
    elseif i == 391
      return rand(p.l391:p.u391)
    else 
      return 0
    end 
  end 


  function PatExp(priv::patientrange, med::patientrange)
    outp_p::patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    outp_m::patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
      # private expected output
    outp_p.count385 = (priv.u385-priv.l385)/2
    outp_p.count386 = (priv.u386-priv.l386)/2
    outp_p.count387 = (priv.u387-priv.l387)/2
    outp_p.count388 = (priv.u388-priv.l388)/2
    outp_p.count389 = (priv.u389-priv.l389)/2
    outp_p.count390 = (priv.u390-priv.l390)/2
    outp_p.count391 = (priv.u391-priv.l391)/2
      # medicaid expected output
    outp_m.count385 = (med.u385-med.l385)/2
    outp_m.count386 = (med.u386-med.l386)/2
    outp_m.count387 = (med.u387-med.l387)/2
    outp_m.count388 = (med.u388-med.l388)/2
    outp_m.count389 = (med.u389-med.l389)/2
    outp_m.count390 = (med.u390-med.l390)/2
    outp_m.count391 = (med.u391-med.l391)/2
    return outp_p, outp_m
  end 






  # 04 02 17 - was "type."  
  immutable coefficients
    distance::Float64
    distsq::Float64
    inten::Float64
    inter::Float64
    distbed::Float64
    closest::Float64
    # can add extras
  end


  type zipcode{T<:Fac}
   code::Int64
   phr::Int64 # may have coefficients differing by PHR
   facilities::Dict{Int64, T}
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
   zips::Dict{Int64, zipcode}
  end

  # Counterfactual-related items:

  type hyrec # quantities of interest within a hospital-year.
    fid::Int64 # hosp ID.
    totbr::Array{Int64,1}  #births
    totlbw::Array{Int64,1}  #lbw births
    totvlbw::Array{Int64,1}  #vlbw births
    deaths::Array{Float64,1}   # deaths
    profit::Array{Float64,1} # hospital revenue.
  end

  type simrun
    # contains the results of some T period sim.
    fips::Int64
    hosprecord::Dict{Int64, hyrec} # track patient volumes and deaths.
    yeartot::Float64
    hasfac::Int64
  end

#TODO - think about changing values to Dict{Array{Int64,1}, simrun} - this will more easily account for the assignment of pairs to have high level facs.

  type mkthistory
    fips::Int64
    # for each FID in the fips, one of these
    values::Dict{Int64, simrun} # note this change.
  end

  type counterhistory
    # One of these for each fips code
    hist::Dict{Int64, mkthistory}
  end


"""
`type nlrec`
- aw::Dict{Int64,Float64}
- psi::Array{Float64,2}
- counter::Dict{Int64, Int64}
First element stores {Action, Continuation Value Approximations}.  Second stores probabilities.  Third is an {Action, Hits} counter.
Stored in a dictionary under a (neighbors, level) tuple-type key.
"""
type nlrec
  aw::Dict{Int64,Float64}
  psi::Array{Float64,2}
  counter::Dict{Int64, Int64}
end

"""
`type hitcount`
- conf::neighbors
- visits::Dict{Int64, Int64}
"""
type hitcount # this type will record visits to a state-action pair
  conf::neighbors
  visits::Dict{Int64,Int64}
end


"""
`vrecord`
- visited::Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, 1} # neighbors, level, action
- totalcnt::Int64 # count all states
Records which state-action pairs have been visited, plus the total iteration count.  The array length is fixed at
1 million so only those state-action pairs are checked.
"""
type vrecord
  visited::Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, 1} # neighbors, level, action
  totalcnt::Int64 # count all states
end


"""
`allvisits`
- all::Dict{Int64, vrecord}
Dictionary of {Fid, VRecord}, which records visits to state-action pairs.
"""
type allvisits
  all::Dict{Int64, vrecord}
end



"""
`type history`
-  path::Dict{neighbors, hitcount}
-  totalcount::Int64
"""
type history
  path::Dict{neighbors, hitcount}
  totalcount::Int64 # records total number of iterations
end


"""
`type shortrec<:ProjectModule.Fac`
-  fid::Int64
-  lat::Float64
-  long::Float64
-  level::Int64
-  truelevel::Int64
-  beds::Int64
-  ns::neighbors
-  choices::Array{Int64, 2}
-  chprobs::WeightVec
-  tbu::Bool
"""
type shortrec<:ProjectModule.Fac
  # needs to take location.  must measure distances.  need to update all of these every time, probably.
  fid::Int64
  lat::Float64
  long::Float64
  level::Int64
  truelevel::Int64
  beds::Int64
  ns::neighbors
  choices::Array{Int64, 2}
  chprobs::WeightVec
  tbu::Bool
end

"""
`type cpats`
-  zp::Int64
-  lat::Float64
-  long::Float64
-  putils::Array{Float64,2}
-  mutils::Array{Float64,2}
-  pwtp::Array{Float64,2}
-  facs::Array{shortrec,1}
-  pcounts::patientcount
-  mcounts::patientcount
"""
type cpats
  zp::Int64
  lat::Float64
  long::Float64
  putils::Array{Float64,2}
  mutils::Array{Float64,2}
  pwtp::Array{Float64,2}
  facs::Array{shortrec,1}
  pcounts::patientcount
  mcounts::patientcount
end


"""
`type cmkt`
-  fid::Int64
-  m::Array{cpats,1}
"""
type cmkt
  fid::Int64
  m::Array{cpats,1}
end


#NB: When the firm exits, can probably restart from the beginning, but keeping the elements in "visited".  We can keep approximating them.
"""
`type simh<:ProjectModule.Fac`
-  fid::Int64
-  lat::Float64
-  long::Float64
-  level::Int64
-  previous::Int64
-  actual::Int64
-  beds::Int64
-  cns::neighbors # must know what current neighbors look like.
-  visited::Dict{nl, nlrec} #possible given "isequal" and "hash" extended for "neighbors"
-  ns::Array{shortrec, 1}
-  mk::cmkt # putting the cmkt into the simh record itself.
-  exit::Bool # did it exit?
-  tbu::Bool # Does the record need to be updated ?
-  converged::Bool # has the hospital converged or not?
"""
type simh<:ProjectModule.Fac
  fid::Int64
  lat::Float64
  long::Float64
  level::Int64
  previous::Int64
  actual::Int64
  beds::Int64
  cns::neighbors # must know what current neighbors look like.
  nfids::Array{Int64,1} # add neighbor fids.  
  visited::Dict{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, nlrec} #change key to tuple of int64 composed of neighbors and level.
  ns::Array{shortrec, 1}
  mk::cmkt # putting the cmkt into the simh record itself.
  exit::Bool
  tbu::Bool
  converged::Bool
end

"""
`type DynState` # this should hold a collection of ALL of the hospitals, for the approximation.
-  all::Array{simh, 1}
"""
type DynState # this should hold a collection of ALL of the hospitals, for the approximation.
  all::Array{simh, 1}
end




  include("LogitEst.jl")
  include("Distance.jl")
  include("DataStructs.jl")
  include("Counter1.jl")
  include("Counter2.jl")
  #include("Reboot.jl")



  export FindUndone
  export CMakeIt
  export SetLevel
  export FillState
  export PatientDraw
  export AllMortality
  export VolMortality
  export LogitEst
  export logitest
  export states1
  export states2
  export states3
  export poly
  export NewPoly
  export NeighborsProd
  export BasicProd
  export perturb
  export distance
  export MakeIt
  export TXSetup
  export ExpandDict
  export MakeNew
  export CreateEmpty
  export MarketPrint
  export NeighborsPrint
  export FacPrint
  export NewEntrantLocation
  export MktSize
  export ChoicesAvailable
  export LevelFunction
  export NeighborAppend
  export NeighborRemove
  export NeighborClean
  export NeighborFix
  export StrictCountyNeighborFix
  export HospFindFirst
  export FidFindFirst
  export MarketCleaner
  export HospUpdate
  export HospPerturb
  export CreateZips
  export FillPPatients
  export FillMPatients
  export FillPatients
  export NewPatients
  export PrintZip
  export ComputeDetUtil
  export WhichZips
  export CalcWTP
  export WTPMap
  export WriteWTP
  export UpdateDeterministic
  export GenPChoices
  export GenMChoices
  export PHistoryAdd
  export MHistoryAdd
  export PDemandMap
  export MDemandMap
  export HospitalClean
  export Restore
  export NewSim
  export Termination
  export PSim
  export TransitionGen
  export CondSum
  export DemandCheck
  export ResultsOut
  export OuterSim
  export AddOO
  export FindFids
  export InitChoice
  export NewPatientsTest
  export DicttoVec
  export PatientCountOut
  export UMap
  export DV
  export NewHospDict
  export PatientsClean
  export ChoiceVector
  export UseThreads
  export ResVec
  export CleanLevelHistory
  export CleanProbHistory
  export CleanWTPHistory
  export CleanDemandHistory
  export RemoveEntrant
  export RecordCopy
  export DrawPatients 
  export PatExp

  # Export Counter1.jl Functions 
  export CategoryReminder
  export Payoff 
  export TermFl
  export FindVLBW
  export MeanCost 
  export CounterSim 
  export Baseline 
  export BaselineCheck
  export DeathCheck 
  export BestOutcome 
  export SimpleResultsPrint 
  export DemandChangeATX 
  export RunCounter
  export CounterClean 
  export CounterCleanResults
  export HHI 
  export RunCounter1
  export ProfitChange
  export ProfitMeanVar


  # Export Counter2.jl functions 
  export KeyCreate 
  export CounterObjects
  export DynStateCreate
  export DynPatients 
  export DetUtils 
  export CounterWTP
  export DSim 
  export DS2 
  export UpdateDUtil
  export HUtil 
  export UpdateCheck 
  export NCheck 
  export GetProbs 
  export FixMainNs
  export FixNN 
  export PatientFind 
  export GetProbCheck 
  export FindWTP 
  export CatchWTP 
  export CatchWTPAll
  export StartingVals
  export ProbUpdate 
  export PolicyUpdate
  export WeightedProbUpdate 
  export MD 
  export DA 
  export WProb 
  export ContError 
  export WhyNaN
  export SinglePay
  export ComputeR 
  export RTuple 
  export ValApprox 
  export KeytoTuple 
  export CheckConvergence 
  export LogitCheck 
  export AllDems 
  export Halt 
  export GetChunk 
  export ChooseAction 
  export ExCheck 
  export PrintVisited 
  export DemandGroup 
  export EasyDemand 
  export ThreeDemands 
  export SimpleDemand 
  export AvgAggState 
  export AllAgg 
  export FixNS 
  export NeighborsCheck 
  export VisitPrint 
  export ProbChange 
  export FindScale 
  export ArrayZero
  export WTPNew
  export DSimNew
  export DemComp
  export GetProb 
  export PatientZero
  export TupleSmash
  export EnumerLevel
  export StateEnumerate
  export DictCopy
  export DictClean
  export ExactConvergence
  export ExactVal
  export EnumUp
  export HospitalDemand
  export StateKey
  export FindComps
  export PatientRev 
  export UpdateD 
  export UtilUp
  export UtilDown 
  export ExactChoice
  export ContVal
  export ContProbs
  export StateRecord 
  export TAddLevel 
  export CombineV
  export PMatch
  export TotalCombine
  export CombineVInput
  export GenStates
  export PartialCombine 
  export TupletoCNS
  export CompsDict
  export NStateKey
  export ConvTest
   
  


  # Export types
  export EntireState
  export chospital
  export LBW
  export patientcollection
  export zipcode
  export coefficients
  export patientcount
  export EntireState
  export Market
  export hospital
  export neighbors
  export DemandHistory
  export initial
  export WTP
  export counterhistory
  export mkthistory
  export simrun
  export hyrec
  export patientrange
    # Counter2.jl types 
  export nlrec
  export hitcount 
  export allvisits
  export history 
  export shortrec 
  export cpats 
  export cmkt 
  export simh 
  export DynState 
  export vrecord 
  export hitcount
  export nlrec

  println("Loaded Module")


  include("Reboot.jl")



end #end of module





####
