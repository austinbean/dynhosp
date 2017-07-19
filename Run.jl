

#addprocs(30) # will work on Amazon.
#addprocs() # will work on home machine.
addprocs(68) # will work on Stampede

loc = pwd()

using ProjectModule

ds1 = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
println(ds1)

res1, res2 = CombinedSim(68; T1 = 40)

println("Estimation Done - Saving")

ds = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
println(ds)

writecsv(loc*"/longobjective"*ds*".csv", res1)
writecsv(loc*"/shortobjective"*ds*".csv", res2)


println("Done Saving")
