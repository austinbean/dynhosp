
#lis = addprocs(1)

module ProjectModule
  using DataFrames # this needs to be inside and outside the module?
  using Distributions

  dir = pwd()
  global pathdata = "";global pathpeople = "";global  pathprograms = "";
  if dir == "/Users/austinbean/Desktop/dynhosp"
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

  abstract Fac


  type hospital <: Fac
    fid::Int64
    lat::Float64
    long::Float64
    name::AbstractString
    fipscode::Int64
    level::Int64
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


  type coefficients
    distance::Float64
    distsq::Float64
    inten::Float64
    inter::Float64
    distbed::Float64
    closest::Float64
    # can add extras
  end

  type zip{T<:Fac}
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
   zips::Dict{Int64, zip}
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



  include("LogitEst.jl")
  include("Distance.jl")
  include("DataStructs.jl")

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
  export EntireState
  export chospital
  export LBW
  export patientcollection
  export zip
  export coefficients
  export patientcount
  export EntireState
  export Market
  export hospital
  export neighbors
  export DemandHistory
  export WTP
  export counterhistory
  export mkthistory
  export simrun
  export hyrec
  export AddOO
  export FindFids
  export InitChoice
  export NewPatientsTest
  export DicttoVec
  export PatientCountOut
  println("Loaded Module")


  include("Reboot.jl")



end #end of module





####
