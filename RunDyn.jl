# Runs the dynamic computation
# Note the parameters wallh and wallm in ExactControl - these determine how long the simulation should take.  

addprocs(68) # will work on Stampede

#addprocs(2) # for testing at home.


loc = pwd()
println("current location: ", loc)
using ProjectModule

facs = ArgVec(ARGS) # take command line input - a list of fids.  

ds1 = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
println("Sim Starting at: ", ds1)

dyn = CounterObjects(1);
res1 = Dict{Int64,Dict{NTuple{10,Int64},Float64}}()
ExactControl(dyn, 11, 30, 1_000_000, 1_000_000, facs; results = res1) # argument is: dynstate, wallhours, wallminutes, exactitlim, approxitlim, list of facs


outp = ResultsWrite(res1)


ds = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")

println("Finished at: ", ds)

writecsv(loc*"/dynresults"*ds*".csv", outp)

println("Done Saving")

