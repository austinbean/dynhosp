

#addprocs()

loc = pwd()

using ProjectModule

res = OuterSim(500; T1 = 50)

res2 = convert(DataFrames.DataFrame, res)

println("Estimation Done - Saving")

DataFrames.writetable(loc*"simulationresults.csv", res2)

println("Done Saving")
