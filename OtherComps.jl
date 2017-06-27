#Other paper computations: 

#=

Miscellaneous model computations for the paper:

=#



"""
`WTPchange(d::DynState)`
Take each firm.  Update the level.  Update the WTP.  Compute the change.

dyn = CounterObjects(50);


"""
function WTPchange(d::DynState)
    outp::Array{Float64,2} = zeros(286, 3)
    for ix in 1:size(dyn.all,1)
        outp[ix,1] = dyn.all[ix].fid
        if dyn.all[ix].level == 1
            println(dyn.all[ix].fid)
            outp[ix,2] = FindWTP(dyn.all[ix])
            println(FindWTP(dyn.all[ix]))
            dyn.all[ix].level = 3
            dyn.all[ix].tbu = true 
            # Ok - this doesn't update yet.  
            UpdateDUtil(dyn.all[ix]) 
            outp[ix,3] = FindWTP(dyn.all[ix])
            println(FindWTP(dyn.all[ix]))
            dyn.all[ix].level = 1
            dyn.all[ix].tbu = true 
            UpdateDUtil(dyn.all[ix])
            println("********")
        elseif dyn.all[ix].level == 2
            outp[ix,2] = FindWTP(dyn.all[ix])
            dyn.all[ix].level = 3 
            dyn.all[ix].tbu = true 
            UpdateDUtil(dyn.all[ix])
            outp[ix,3] = FindWTP(dyn.all[ix])
            dyn.all[ix].level = 2
            dyn.all[ix].tbu = true 
            UpdateDUtil(dyn.all[ix])
        elseif dyn.all[ix].level == 3
            # do nothing.
        end 
    end 
    return outp
end 

"""
`DUtilCheck`
See whether the utilities are updated correctly for "shortrecs"
"""
function DUtilCheck(h::simh)
    println("Facilities")
    for el in h.mk.m 
        println(el.zp, " ", size(unique(el.putils[1,:]))[1], "  ", size(unique(el.mutils[1,:]))[1])
    end 


    println("Own Utility Update Check: ")
    for el in h.mk.m 
        println(el.zp, "  ", el.putils[2, findin(el.putils[1,:], h.fid)])
    end 
    h.level = 3
    h.tbu = true 
    UpdateDUtil(h)
    println("********")
    for el in h.mk.m 
        println(el.zp, "  ", el.putils[2, findin(el.putils[1,:], h.fid)])
    end
    h.level = 1 
    h.tbu = true 
    UpdateDUtil(h) 
end 

#=
This computation relates to determining the market sizes of individual hositals by checking the total number of births in 
all zip codes attached to that hospital every year.  It is merged with some data from Stata.  Specifically, see TX Patient Uncertainty.do

=#
mx = 0;
for el in dyn.all 
    if size(el.mk.m,1) > mx 
        mx = size(el.mk.m,1)
    end 
end 
println("Max zips = ", mx)

zps = zeros(Int64, mx+1, size(dyn.all,1) ); # NB - each COLUMN is an fid 

for el in 1:size(dyn.all,1) # COLUMNS
    zps[1, el] = dyn.all[el].fid  
    for z in 1:size(dyn.all[el].mk.m,1) # ROWS 
        zps[z+1,el] = dyn.all[el].mk.m[z].zp # zip 
    end 
end 
zps = transpose(zps)
writecsv( "/Users/austinbean/Desktop/dynhosp/hospzips.csv", zps)




#=

typetest 


=#

mutable struct WT 
    w::Array{Float64,1}
end 

struct WTPee
    w::Array{Float64, 1}
end

w1 = WTPee(Array{Float64,1}(50))
w2 = WT(Array{Float64,1}(50))



for i = 1:10
@time for k1 = 1:50 
    w1.w[k1] = 3.0
end 
for k = 1:50 
    w1.w[k] = 0.0
end 
end 

for i = 1:10
@time for k1 = 1:50 
    w2.w[k1] = 3.0
end 
for k = 1:50 
    w2.w[k] = 0.0
end 
end 



struct w1
    d::Dict{Int64,Int64}
end 

ab1 = WTPee(Array{Float64,1}())

ab2 = WTPee(Array{Float64,1}(50))

@time for el =1:50
    push!(ab1.w, rand())
end 

@time for el =1:50 
    ab2.w[el] = rand()
end 

for j = 1:10
@time for i = 1:50
 ab2.w[i] = i 
end 
for k = 1:50
    ab2.w[k] = 0.0
end 
end 


ab = w1(Dict{Int64,Int64}())


for el in Texas.ms
    for k in el.config 
        levl = (-1,-1)
         if k.level == 1
           levl = (0,0)
         elseif k.level == 2
           levl = (1,0)
         elseif k.level == 3
           levl = (0,1)
         end
        n1, n2, n3 = MktSize(k.neigh)
        nv = [k.neigh.level105; k.neigh.level205; k.neigh.level305; k.neigh.level1515; k.neigh.level2515; k.neigh.level3515; k.neigh.level11525; k.neigh.level21525; k.neigh.level31525 ]
        try a1 = logitest(levl, n1, n2, n3, nv)
        catch er1
            println(er1, "the error")
            if isa(er1,ProjectModule.ValueException)
                println("didn't work:")
                println(levl)
                println(n1, " ", n2, " ", n3)
                println(nv)
                println(typeof(nv))
            else 
                # do nothing
            end 
        end 
    end 
end 



#=
import Base.+
Base.+(a1::String, a2::String) = string(a1,a2)
=#





# TODO - must get the output of ResultsOut and ResultsOutVariant into the @parallel.  

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);

outp1, outp2 = @sync @parallel (+) for j = 1:MCcount
    # probably need a new function to make this work.  

end 


# Some Benchmarks related to PSim()

@benchmark WriteWTP($WTPMap(patients, Texas), $Texas, 1)
BenchmarkTools.Trial:
  memory estimate:  2.56 MiB
  allocs estimate:  14741
  --------------
  minimum time:     6.164 ms (0.00% GC)
  median time:      7.412 ms (0.00% GC)
  mean time:        7.681 ms (3.55% GC)
  maximum time:     12.744 ms (31.44% GC)
  --------------
  samples:          650
  evals/sample:     1

@benchmark WTPMap($patients, $Texas)
BenchmarkTools.Trial:
  memory estimate:  2.50 MiB
  allocs estimate:  10737
  --------------
  minimum time:     5.307 ms (0.00% GC)
  median time:      5.828 ms (0.00% GC)
  mean time:        6.185 ms (3.11% GC)
  maximum time:     10.855 ms (26.12% GC)
  --------------
  samples:          807
  evals/sample:     1

@benchmark NewHospDict($Texas)
BenchmarkTools.Trial:
  memory estimate:  42.53 KiB
  allocs estimate:  304
  --------------
  minimum time:     18.695 μs (0.00% GC)
  median time:      21.220 μs (0.00% GC)
  mean time:        36.018 μs (14.99% GC)
  maximum time:     5.607 ms (97.76% GC)
  --------------
  samples:          10000
  evals/sample:     1

@benchmark GenPChoices($patients, $d1, $arry1)
BenchmarkTools.Trial:
  memory estimate:  33.44 MiB
  allocs estimate:  644479
  --------------
  minimum time:     143.279 ms (1.80% GC)
  median time:      147.961 ms (3.44% GC)
  mean time:        148.439 ms (2.72% GC)
  maximum time:     153.695 ms (3.65% GC)
  --------------
  samples:          34
  evals/sample:     1


#=

=#