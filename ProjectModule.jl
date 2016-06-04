# Module for Main Work:
# push this change to the LOAD_PATH: push!(LOAD_PATH, "/Users/austinbean/Desktop/dynhosp")
#
# using DataFrames
# @everywhere using DataFrames
#
# using Distributions
# @everywhere using Distributions


module ProjectModule
  using DataFrames # this needs to be inside and outside the module?
  using Distributions 
  include("/Users/austinbean/Desktop/dynhosp/LogitEst.jl")
  include("/Users/austinbean/Desktop/dynhosp/Distance.jl")
  include("/Users/austinbean/Desktop/dynhosp/Simulator.jl")
  include("/Users/austinbean/Desktop/dynhosp/PerturbSimulation.jl")
  include("/Users/austinbean/Desktop/dynhosp/DynamicValue.jl")
  include("/Users/austinbean/Desktop/dynhosp/DemandModel.jl")
  include("/Users/austinbean/Desktop/dynhosp/MainFunction.jl")
  include("/Users/austinbean/Desktop/dynhosp/PerturbAction.jl")
  include("/Users/austinbean/Desktop/dynhosp/ParallelFunctions.jl")

  export LogitEst, Simulator, PerturbSimulator, logitest, states1, states2, states3, poly, Mainfun, DemandModel, EntrantsU, rowchange, rowfindfid, fidfinder, sendto, passobj, perturb, distance, DynamicValue

end #end of module


#=


# Now send to all workers:
sendto(lis; data=data)



# Send to all workers:
sendto(lis; people=people)

sendto(lis, demandmodelparameters = demandmodelparameters)

=#
