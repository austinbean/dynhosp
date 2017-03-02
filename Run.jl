

#addprocs(30) # will work on Amazon.
#addprocs() # will work on home machine.

loc = pwd()

using ProjectModule

res = OuterSim(6; T1 = 6)

println("Estimation Done - Saving")

ds = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")

writecsv(loc*"/simulationresults"*ds*".csv", res)

println("Done Saving")
