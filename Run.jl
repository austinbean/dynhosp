

#addprocs(30) # will work on Amazon.
#addprocs() # will work on home machine.
#addprocs(68) # will work on Stampede
loc = pwd()
println("current location: ", loc)
using ProjectModule

ds1 = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
println("Sim Starting at: ", ds1)

res1, res2 = CombinedSim(1; T1 = 1) # 476/68 = 7 reps/core.  

println("Estimation Done - Saving")
println("size? ", size(res1), " ", size(res2))

ds = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")

println("Finished at: ", ds)

println("locs: ", loc*"/longobjective"*ds*".csv")
println("locs2: ", loc*"/shortobjective"*ds*".csv")

#writecsv(loc*"/longobjective"*ds*".csv", res1)
#writecsv(loc*"/shortobjective"*ds*".csv", res2)


println("Done Saving")
