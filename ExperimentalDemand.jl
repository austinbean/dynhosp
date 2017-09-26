
"""

`ChoiceVectorEXP(pd::Dict{Int64, Float64},
                         dt::Dict{Int64, patientcount},
                         ch::Array{Int64,1},
                         x::patientcount,
                         cts::Dict{Int64,Int64})`

STILL AN EXPERIMENT


Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
inpt = ones(Int64, 1550); # largest group is 1511
fids1, utils1 = DV(patients.zips[78759].pdetutils)
tar = zeros(utils1)
# nm is patients.zips[78759].ppatients.count391
UseThreads(inpt, fids1, utils1,  500)
TODO - fix this - dic1 does not exist yet.
newd = Dict{Int64,Int64}()
for k1 in keys(dic1)
  newd[k1] = 0
end 

@benchmark ChoiceVectorEXP(patients.zips[78759].pdetutils, dic1, inpt, patients.zips[78759].ppatients, newd)

BenchmarkTools.Trial:
  memory estimate:  608 bytes
  allocs estimate:  10
  --------------
  minimum time:     581.567 μs (0.00% GC)
  median time:      585.169 μs (0.00% GC)
  mean time:        603.198 μs (0.00% GC)
  maximum time:     1.334 ms (0.00% GC)
  --------------
  samples:          8273
  evals/sample:     1

"""
function ChoiceVectorEXP(pd::Dict{Int64, Float64},
                         dt::Dict{Int64, patientcount},
                         ch::Array{Int64,1},
                         x::patientcount,
                         cts::Dict{Int64,Int64})
  #fids::Array{Int64,1}, utils::Array{Float64,1} = DV(pd) # very quick ≈ 300 ns -  this counts for 4 allocations.  
  # NB - the above is unnecessary with the revised function UMapDict.  
  loc::Int64 = 1
  for nm in x
    UseThreads(ch, utils, nm)                      # ≈ 179 μs, for nm = 300, ≈ 12.504 μs for nm = 20
    Frequency(cts, ch, nm)                               # compute the frequencies of the choices.
    if loc == 1
      for i in keys(cts) 
        dt[i].count385 = cts[i]                          
      end
      loc += 1
    elseif loc==2
      for i in keys(cts)
        dt[i].count386 = cts[i]
      end
      loc += 1
    elseif loc==3
      for i in keys(cts)
        dt[i].count387 = cts[i]
      end
      loc += 1
    elseif loc==4
      for i in keys(cts)
        dt[i].count388 = cts[i]
      end
      loc += 1
    elseif loc==5
      for i in keys(cts)
        dt[i].count389 = cts[i]
      end
      loc += 1
    elseif loc==6
      for i in keys(cts)
        dt[i].count390 = cts[i]
      end
      loc += 1
    elseif loc==7
      for i in keys(cts)
        dt[i].count391 = cts[i]
      end
      loc += 1
    else 
      # do nothing.  
    end
    ResVec(ch) #reset the vector. memory estimate:  0 bytes - 471.337 ns (0.00% GC)
  end
end





"""
`UseThreadsEXP()`

THIS IS AN EXPERIMENT.

This would take the dict directly, instead of splitting it into two vectors first.  

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
inpt = ones(Int64, 1550); # largest group is 1511

UseThreadsEXP(inpt, patients.zips[78759].pdetutils, 500)

BenchmarkTools.Trial:
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     442.187 μs (0.00% GC)
  median time:      445.259 μs (0.00% GC)
  mean time:        462.430 μs (0.00% GC)
  maximum time:     1.110 ms (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1

"""
function UseThreadsEXP(inpt::Array{Int64,1},utils::Dict{Int64,Float64}, x::Int64)
  #Threads.@threads # using this DOES allocate.  If it isn't faster, cut it out.  
  # TODO - can this use the dict instead of fids, utils arrays?  This would take EITHER changes to UMAP or 
  # just call this with allocated vectors...?
  for i = 1:x
    inpt[i] = UMapDict(utils)
  end
end


"""
`Frequency`

STILL AN EXPERIMENT.

This counts the frequency of all elements in the vector arr of Ints.
It requires a dict of Int64's which will be {fid, count}
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);
inpt = ones(Int64, 1550); # largest group is 1511
fids1, utils1 = DV(patients.zips[78759].pdetutils)
tar = zeros(utils1)
# nm is patients.zips[78759].ppatients.count391
UseThreads(inpt, fids1, utils1,  500)


newd = Dict{Int64,Int64}()
for k1 in keys(dic1)
  newd[k1] = 0
end 

@benchmark Frequency(newd,  inpt, 500)
BenchmarkTools.Trial:
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     28.626 μs (0.00% GC)
  median time:      28.802 μs (0.00% GC)
  mean time:        28.832 μs (0.00% GC)
  maximum time:     74.631 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
"""
function Frequency(pds::Dict{Int64,Int64}, arr::Array{Int64,1}, nm::Int64)
  for el in keys(pds)        # this will operate with a smaller set of keys 
    pds[el] = 0
  end 
  for i = 1:nm               # only go over elements which are used in the vector arr.  
    if arr[i] != 1           # skip the default in the vector 
      if haskey(pds, arr[i]) # this will always have the keys 
        pds[arr[i]] += 1
      else 
        pds[arr[i]] = 1
      end 
    end 
  end 
end 


"""
`UMapDict`
Experimental version of umap taking a dict.

This is about half the speed, but does not allocate.
The advantage of replacing this in the function ChoiceMap is that the step 
DV = ..., ... which allocates two arrays could be avoided.  This is about half the 
allocations in that function.  

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = NewPatients(Texas);

UMapDict(patients.zips[78702].pdetutils)

BenchmarkTools.Trial:
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     904.737 ns (0.00% GC)
  median time:      925.566 ns (0.00% GC)
  mean time:        961.796 ns (0.00% GC)
  maximum time:     2.771 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     38

"""
function UMapDict(utils::Dict{Int64,Float64})::Int64
  const dist_μ::Int64 = 0
  const dist_σ::Int64 = 1
  const dist_ξ::Int64 = 0
  d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
  maxm::Int64 = 1                   # this is an index to a fid.  Will record max index.
  maxu::Float64 = 0.0               # holds the value of max utility 
  tem::Float64 = 0.0                # holds interim utility value
  for i in keys(utils)
    tem = utils[i]+rand(d)
    if tem>maxu  
        maxm = i                    # replace index if greater
        maxu = tem                  # replace max util 
    end 
  end 
  return maxm::Int64          # return fid of max util value. 
end 
