


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
  end


  include(pathprograms*"LogitEst.jl")
  include(pathprograms*"Distance.jl")
  #include(pathprograms*"Simulator.jl")
  #include(pathprograms*"PerturbSimulation.jl")
  #include(pathprograms*"DynamicValue.jl")
  #include(pathprograms*"DemandModel.jl")
  #include(pathprograms*"MainFunction.jl")
  #include(pathprograms*"PerturbAction.jl")
  #include(pathprograms*"ParallelFunctions.jl")
  #include(pathprograms*"ParMainFun.jl")
  #include(pathprograms*"WTP.jl")



  export LogitEst, logitest, states1, states2, states3, poly, perturb, distance
  println("Loaded Module")
end #end of module





####
