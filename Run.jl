

#addprocs()

loc = pwd()

using ProjectModule

res = OuterSim(4; T1 = 1)

res2 = convert(DataFrame, res)

writetable(loc*"simulationresults.csv", res2)
