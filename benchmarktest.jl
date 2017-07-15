importall ProjectModule

using BenchmarkTools 

println("Data Objects Setup: ")
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
dic1 = NewHospDict(Texas);
inpt = ones(Int64, 1550); # largest group is 1511
d1 = WTPMap(patients, Texas);
fids1, utils1 = DV(patients.zips[78759].pdetutils);
ut = [0.1, 0.2, 0.3, 0.4, 0.5];
fi = [111, 222, 333, 444, 19];

println("Number of Threads: ", Threads.nthreads() )

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
    println("memory estimate: ", x.memory)
    println("allocations estimate: ", x.allocs)
    println("--------")
    println("minimum time: ", minimum(x.times))
    println("median time: ", median(x.times))
    println("mean time: ", mean(x.times))
    println("maximum time: ", maximum(x.times))
    println("--------")
    println("samples: ", x.params.samples)
    println("        ")
end 


for k1 in keys(res)
    println("Benchmark:         ", k1)
    BRPrint(res[k1])
end 


#=
    ME:
Benchmark:         MDemandMap
memory estimate: 22944
allocations estimate: 573
--------
minimum time: 77055.0
median time: 78440.0
mean time: 89145.5664
maximum time: 8.039676e6
--------
samples: 10000

    STAMPEDE:
Benchmark:         MDemandMap
memory estimate: 22944
allocations estimate: 573
--------
minimum time: 690533.0
median time: 754982.0
mean time: 779351.5574
maximum time: 6.3525857e7
--------
samples: 10000


######################

    ME:    
Benchmark:         UpdateDeterministic
memory estimate: 7492576
allocations estimate: 453242
--------
minimum time: 2.0536423e7
median time: 2.2822403e7
mean time: 2.383370995952381e7
maximum time: 3.9156113e7
--------
samples: 10000

    STAMPEDE:
Benchmark:         UpdateDeterministic
memory estimate: 7492576
allocations estimate: 453242
--------
minimum time: 1.20593127e8
median time: 1.22866134e8
mean time: 1.2802742982278481e8
maximum time: 1.45036322e8
--------
samples: 10000

#######################

    ME:

Benchmark:         UMap
memory estimate: 0
allocations estimate: 0
--------
minimum time: 284.0
median time: 321.0
mean time: 341.4901
maximum time: 53706.0
--------
samples: 10000

    STAMPEDE:

Benchmark:         UMap
memory estimate: 0
allocations estimate: 0
--------
minimum time: 1704.0
median time: 1794.0
mean time: 1873.0699
maximum time: 41667.0
--------
samples: 10000

########################
    ME:
Benchmark:         UseThreads
memory estimate: 48
allocations estimate: 1
--------
minimum time: 177540.0
median time: 179971.0
mean time: 190044.4739
maximum time: 1.967888e6
--------
samples: 10000

    STAMPEDE:
Benchmark:         UseThreads
memory estimate: 48
allocations estimate: 1
--------
minimum time: 1.052319e6
median time: 1.059074e6
mean time: 1.0790485258414766e6
maximum time: 1.1073888e7
--------
samples: 10000


#########################

    ME: 

Benchmark:         GenPChoices
memory estimate: 6869379
allocations estimate: 242708
--------
minimum time: 1.35636428e8
median time: 1.46256591e8
mean time: 1.5215130486363637e8
maximum time: 2.66270969e8
--------
samples: 10000

    STAMPEDE:

Benchmark:         GenPChoices
memory estimate: 6865052
allocations estimate: 242433
--------
minimum time: 7.96608235e8
median time: 8.015005365e8
mean time: 8.780128120833334e8
maximum time: 1.632352646e9
--------
samples: 10000

##########################

    ME:
Benchmark:         ChoiceVector
memory estimate: 19728
allocations estimate: 1170
--------
minimum time: 301343.0
median time: 305209.5
mean time: 333479.9599
maximum time: 8.420701e6
--------
samples: 10000

    STAMPEDE:

Benchmark:         ChoiceVector
memory estimate: 18624
allocations estimate: 1101
--------
minimum time: 1.732087e6
median time: 1.763049e6
mean time: 1.8018968434955673e6
maximum time: 1.9188356e7
--------
samples: 10000

################################

    ME:

Benchmark:         GenMChoices
memory estimate: 11590576
allocations estimate: 607959
--------
minimum time: 1.53599632e8
median time: 1.64535552e8
mean time: 1.7081858683050847e8
maximum time: 2.42647234e8
--------
samples: 10000

    STAMPEDE:

Benchmark:         GenMChoices
memory estimate: 10831328
allocations estimate: 560506
--------
minimum time: 9.0403674e8
median time: 9.12362025e8
mean time: 9.175241630909091e8
maximum time: 9.29613851e8
--------
samples: 10000

##############################

    ME: 

Benchmark:         WTPMap
memory estimate: 2618048
allocations estimate: 10737
--------
minimum time: 5.304324e6
median time: 5.4884785e6
mean time: 5.867183841950647e6
maximum time: 1.7335866e7
--------
samples: 10000

    STAMPEDE:

Benchmark:         WTPMap
memory estimate: 2618048
allocations estimate: 10737
--------
minimum time: 2.2217959e7
median time: 2.2884036e7
mean time: 2.3363071742990654e7
maximum time: 3.3663264e7
--------
samples: 10000

##########################

    ME:

Benchmark:         PDemandMap
memory estimate: 22944
allocations estimate: 573
--------
minimum time: 76896.0
median time: 79103.5
mean time: 91795.5511
maximum time: 8.056273e6
--------
samples: 10000

    STAMPEDE:

Benchmark:         PDemandMap
memory estimate: 22944
allocations estimate: 573
--------
minimum time: 681052.0
median time: 739877.0
mean time: 765051.9478
maximum time: 6.471412e7
--------
samples: 10000

############################

    ME:

Benchmark:         WriteWTP
memory estimate: 64064
allocations estimate: 4004
--------
minimum time: 343310.0
median time: 350168.5
mean time: 432527.8795
maximum time: 8.264941e6
--------
samples: 10000

    STAMPEDE:

Benchmark:         WriteWTP
memory estimate: 64064
allocations estimate: 4004
--------
minimum time: 1.994092e6
median time: 2.117076e6
mean time: 2.1765138329292485e6
maximum time: 1.7928729e7
--------
samples: 10000

=#





