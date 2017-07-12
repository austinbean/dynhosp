

#addprocs(30) # will work on Amazon.
#addprocs() # will work on home machine.

loc = pwd()

using ProjectModule


res1, res2 = CombinedSim(4; T1 = 2)

println("Estimation Done - Saving")

ds = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")

writecsv(loc*"/longobjectiveresults"*ds*".csv", res1)
writecsv(loc*"/shortobjectiveresults"*ds*".csv", res2)


println("Done Saving")
