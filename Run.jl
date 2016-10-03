

#addprocs()

loc = pwd()

using ProjectModule

res = OuterSim(5; T1 = 10)

res2 = convert(DataFrames.DataFrame, res)

println("Estimation Done - Saving")

DataFrames.writetable(loc*"simulationresults.csv", res2)

println("Done Saving")
