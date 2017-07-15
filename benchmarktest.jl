importall ProjectModule

using BenchmarkTools 

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
dic1 = NewHospDict(Texas);
inpt = ones(Int64, 1550); # largest group is 1511
d1 = WTPMap(patients, Texas);
ids1, utils1 = DV(patients.zips[78759].pdetutils);
ut = [0.1, 0.2, 0.3, 0.4, 0.5];
fi = [111, 222, 333, 444, 19];

bg = BenchmarkGroup()
bg["UMap"] = @benchmarkable UMap($ut, $fi)
bg["GenPChoices"] = @benchmarkable GenPChoices($patients, $dic1, $inpt)
bg["GenMChoices"] = @benchmarkable GenMChoices($patients, $dic1, $inpt)
bg["ChoiceVector"] = @benchmarkable ChoiceVector($patients.zips[78759].pdetutils, $dic1, $inpt, $patients.zips[78759].ppatients)
bg["WTPMap"] = @benchmarkable WTPMap($patients, $Texas)
bg["WriteWTP"] = @benchmarkable WriteWTP($d1, $Texas, 1)
bg["PDemandMap"] = @benchmarkable PDemandMap($dic1, $Texas, 1)
bg["MDemandMap"] = @benchmarkable MDemandMap($dic1, $Texas, 1)
bg["UseThreads"] = @benchmarkable UseThreads($inpt, $fids1, $utils1, 297)
bg["UpdateDeterministic"] = @benchmarkable UpdateDeterministic($patients)
res = run(bg, verbose = true, seconds = 10)


function BRPrint(x::BenchmarkTools.Trial)
    println("memory: ", x.memory)
    println("allocations: ", x.allocs)
    println("--------")
    println("minimum time: ", minimum(x.times))
    println("median time: ", median(x.times))
    println("mean time: ", mean(x.times))
    println("maximum time: ", maximum(x.times))
    println("--------")
    println("        ")
end 


for k1 in keys(res)
    println("Benchmark:         ", k1)
    BRPrint(res[k1])
end 


#=

# things to test:
# important components of NewSim - these are informative about PSim too.
# WriteWTP, GenPChoices, PDemandMap, GenMChoices, MDemandMap, UpdateDeterministic.
# UMap, ChoiceVector, Mapchoices, WriteWTP 
println("number of threads: ", Threads.nthreads() )

println("Benchmarking UMap")
bu = @benchmark UMap($ut, $fi)

# GenPChoices, GenMChoices.
println("Benchmarking Gen P Choices")
bgp = @benchmark GenPChoices($patients, $dic1, $inpt)
println("Benchmarking Gen M Choices")
bgm =@benchmark GenMChoices($patients, $dic1, $inpt) 

# ChoiceVector:
println("Benchmarking Choice Vector")
bcv = @benchmark ChoiceVector($patients.zips[78759].pdetutils, $dic1, $inpt, $patients.zips[78759].ppatients)

# WTPMap 
println("Benchmarking of WTP Map")
bwtpm = @benchmark WTPMap($patients, $Texas)

# WriteWTP 
println("Benchmarking WriteWTP")
bwwtp = @benchmark WriteWTP($d1, $Texas, 1)

# PDemandMap 
println("Benchmarking PDemandMap")
bpdm = @benchmark PDemandMap($dic1, $Texas, 1)

# MDemandMap 
println("Benchmarking MDemandMap")
bmdm = @benchmark MDemandMap($dic1, $Texas, 1)

# UseThreads.
# nm is patients.zips[78759].ppatients.count391
println("Benchmarking UseThreads ")
but = @benchmark UseThreads($inpt, $fids1, $utils1, 297)

# UpdateDeterministic 
println("Benchmarking UpdateDeterministic ")
bud = @benchmark UpdateDeterministic($patients)

=#





