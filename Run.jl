

#addprocs()

loc = pwd()

using ProjectModule

res = OuterSim(1; T1 = 2)

res2 = convert(DataFrames.DataFrame, res)

println("Estimation Done - Saving")

DataFrames.writetable(loc*"simulationresults.csv", res2)

println("Done Saving")
