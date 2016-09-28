

#addprocs()

loc = pwd()

using ProjectModule

res = OuterSim(4; T1 = 2)

res2 = convert(DataFrames.DataFrame, res)

DataFrames.writetable(loc*"simulationresults.csv", res2)
