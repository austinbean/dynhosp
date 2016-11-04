

#addprocs(30) # will work on Amazon.
addprocs() # will work on home machine.

loc = pwd()

using ProjectModule

res = OuterSim(2; T1 = 1)

res2 = convert(DataFrames.DataFrame, res)

println("Estimation Done - Saving")

DataFrames.writetable(loc*"simulationresults.csv", res2)

println("Done Saving")
