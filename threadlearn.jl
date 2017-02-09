



#=
After importing the project module.
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);
GenPChoices(patients, Texas);
NB: DicttoVec(patients.zips[xxxx].pdetutils)

=#


"""
`function UMap(utils::Array{Float64,1}, fids::Array{Int64,1}, temparr::Array{Float64,1})`

This function takes:
- Array of Utils
- Array of Fids
- Temporary Array

Computes utility + random component, maps out corresponding FID.

testing:
ut = [0.1, 0.2, 0.3, 0.4, 0.5];
fi = [111, 222, 333, 444, 19];
ta = [0.0, 0.0, 0.0, 0.0, 0.0];
UMap(ut, fi, ta)
"""
function UMap(utils::Array{Float64,1},
              fids::Array{Int64,1},
              temparr::Array{Float64,1};
              dist_μ = 0,
              dist_σ = 1,
              dist_ξ = 0,
              d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))::Int64
  return fids[indmax(utils+rand!(d, temparr))]
end



"""
`function DV(d::Dict{Int64, Float64})::Tuple{Array{Int64,1},Array{Float64,1}}`
This is a more efficient version of `DicttoVec`.  Takes a dictionary of {Int64,Float64}
and returns two vectors.


#Testing on Choice Data:
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);

DV(patients.zips[78702].pdetutils)
"""
function DV(d::Dict{Int64, Float64})::Tuple{Array{Int64,1},Array{Float64,1}}
  out1::Array{Int64,1} = zeros(Int64, d.count) #for the keys/FIDs
  out2::Array{Float64,1} = zeros(Float64, d.count) #for the utils
  for (i,k) in enumerate(keys(d))
    out1[i] = k
    out2[i] = d[k]
  end
  return out1, out2
end


"""
`function ChoiceVector(utils::Array{Float64,1}, fids::Array{Float64,1}, x::Int64)`
Allocates an array, then uses threading to allocate calls to `UMap` across different threads.  Output
is a Dict{Float64, Int64} of facilities and patient counts.

#Testing on Choice Data:
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);

ChoiceVector(patients.zips[78702].pdetutils, 50)
#NB - there is a performance issue in that the type of many of these is read as Core.Box, but this is #15276 or something.
"""
function ChoiceVector(pd::Dict{Int64, Float64}, x::Int64)
  fids::Array{Int64,1}, utils::Array{Float64,1} = DV(pd)
  outp::Array{Int64,1} = zeros(Int64, x)
  temparry::Array{Float64, 1} = zeros(fids)
  for el in fids
    dt[el] = 0
  end
  Threads.@threads for i = 1:x
    outp[i] = UMap(utils, fids, temparry)
  end
  dt::Dict{Int64, Int64} = Dict()

  for i = 1:x
    dt[outp[i]] += 1
  end
  return dt::Dict{Int64, Int64}
end



"""
`function DictCombine(outp, arg...)`
This function takes an arbitrary number of dictionaries and adds them to outp.
Dicts must be Dict{Int64, Int64}.

D1 = Dict(1 => 3, 2 => 4, 3 => 16, 4 => 12)
D2 = Dict( 6=> 36, 8=>64 , 10=>100 , 12=> 144, 14=>1414 )
D3 = Dict( 7=> 49, 9=> 81, 11=> 121, 13=> 169, 15=> 225)
"""
function DictCombine(outp, arg...)::Dict{Int64, Int64}
  # This could easily be rewritten to take a patientcount as the second arg to the dict.
  for (i,dct) in enumerate(arg) # arg is enumerated, but that generates an (int, dict) pair.
    for k in keys(dct)
      if haskey(outp, k)
        outp[k] += dct[k]
      else
        outp[k] = dct[k]
      end
    end
  end
  return outp
end



# Each zip returns a dictionary of {fid, demand} for ONE DRG.
# work with that simple case and then make it worse.
function attempt2(pz::patientcollection)
  outp::Dict{Float64,Int64} = Dict{Float64, Int64}()
  for ky in keys(pz.zips)
    DictCombine(outp, ChoiceVector(pz.zips[ky].pdetutils, pz.zips[ky].ppatients.count391))
  end
  return outp
end








###




#=
Some investigations of threading and how it works below.

=#



#=
# Curious problem - w/out the loop, all types are inferred correctly.
# w/ the loop and @threads, all become boxed.

function testtype(x::Int64)
  outp::Array{Int64,1} = zeros(Int64, x)
  Threads.@threads for i = 1:x
    outp[i] = i
  end
end

function testtype2(x::Int64)
  outp::Array{Int64,1} = zeros(Int64, x)
  for i = 1:x
    outp[i] = i
  end
end

testtype(1);
testtype2(1);

@code_warntype testtype(1);
@code_warntype testtype2(1);

=#





function thtest(x1::Int64, x2::Int64, x3::Int64)
  outp::threadtest = threadtest(0,0,0)
  temparr::Array{Float64,1} = [0.0; 0.0; 0.0]
  testarr::Array{Float64, 2} = [1.0 2 3; 0.0 0 0]
  count1::Int64 = 0
  count2::Int64 = 0
  count3::Int64 = 0
  for c1 = 1:x1
    count1 += 1

    mind = testarr[1, indmax(testarr[2,:]+randn!(temparr))]
    if mind == 1
      outp.a += 1
    elseif mind == 2
      outp.b += 1
    elseif mind == 3
      outp.c += 1
    else
      println("wat")
    end
  end
  for c2 = 1:x2
    count2 += 1

    mind = testarr[1, indmax(testarr[2,:]+randn(3))]
    if mind == 1
      outp.a += 1
    elseif mind == 2
      outp.b += 1
    elseif mind == 3
      outp.c += 1
    else
      println("wat")
    end
  end
  for c3 = 1:x3
    count3 += 1

    mind = testarr[1, indmax(testarr[2,:]+randn(3))]
    if mind == 1
      outp.a += 1
    elseif mind == 2
      outp.b += 1
    elseif mind == 3
      outp.c += 1
    else
      println("wat")
    end
  end
  return outp, count1, count2, count3
end

function checkthreads()
  return println(Threads.nthreads())
end


function thtest2(x1::Int64, x2::Int64, x3::Int64)
  outp::threadtest = threadtest(0,0,0)
  temparr::Array{Float64,1} = [0.0; 0.0; 0.0]
  testarr::Array{Float64, 2} = [1.0 2 3; 0.0 0 0]
  count1::Int64 = 0
  count2::Int64 = 0
  count3::Int64 = 0
  Threads.@threads for c1 = 1:x1
    mind = testarr[1, indmax(testarr[2,:]+randn!(temparr))]
    count1 += 1
    if mind == 1
      outp.a += 1
    elseif mind == 2
      outp.b += 1
    elseif mind == 3
      outp.c += 1
    else
      println("wat")
    end
  end
  Threads.@threads for c2 = 1:x2
    mind = testarr[1, indmax(testarr[2,:]+randn(3))]
    count2 += 1
    if mind == 1
      outp.a += 1
    elseif mind == 2
      outp.b += 1
    elseif mind == 3
      outp.c += 1
    else
      println("wat")
    end
  end
  #Threads.@threads
  for c3 = 1:x3
    #=
    Ok - threads don't work here.  What's wrong?  It is something like: it's not waiting for all of them to finish.
    OR it doesn't have proper access to this "mind" variable.  Probably this is not someting where access to the
    shared resource is possible?
    Maybe I need as many copies of mind as there are threads?  Try that.
    Options - allocated according to number of threads or allocate according to number of iterations x3?
    =#
    mind = testarr[1, indmax(testarr[2,:]+randn(3))]
    count3 += 1
    if mind == 1
      outp.a += 1
    elseif mind == 2
      outp.b += 1
    elseif mind == 3
      outp.c += 1
    else
      println("wat")
    end
  end
  return outp, count1, count2, count3
end


function countit(arr::Array{Int64,1})
  outp::threadtest = threadtest(0, 0, 0)
  for el in arr
    if el == 1
      outp.a += 1
    elseif el == 2
      outp.b += 1
    elseif el == 3
      outp.c += 1
    end
  end
  return outp
end

function V1(x::Int64)
#  outp::threadtest = threadtest(0,0,0)
  temparr::Array{Float64,2} = [0.0 0.0 0.0]
  testarr::Array{Float64, 2} = [0.0 0.0 0.0]
  a = zeros(Int64, x)
#  t2::Array{Float64, 2} = repmat(testarr[2,:], x, 1)   # only need to replicate second row
  for k1 = 1:x
    a[k1] = mapper(testarr, temparr)
  end
  return countit(a)
end


function V2(x::Int64)
#  outp::threadtest = threadtest(0,0,0)
  temparr::Array{Float64,2} = [0.0 0.0 0.0]
  testarr::Array{Float64, 2} = [0.0 0.0 0.0]
  a = zeros(Int64, x)
#  t2::Array{Float64, 2} = repmat(testarr[2,:], x, 1)   # only need to replicate second row
  Threads.@threads for k1 = 1:x
    a[k1] = mapper(testarr, temparr)
  end
  return countit(a)
end

function F1(x::Int64)
  outp = zeros(Int64, x)
  for k1 = 1:x
    outp[k1] = k1^2
  end
  return outp
end

function F2(x::Int64)
  outp = zeros(Int64, x)
  Threads.@threads for k1 = 1:x
    outp[k1] = k1^2
  end
  return outp
end

#=
- The key should be to make sure that the ORDER is preserved.
- But once the order is preserved, the utilities work in the way expected and
this is pretty quick, I think.

=#


function doit()
  println("Number of threads")
  checkthreads()
  println("testing ")
  thtest(1, 1, 1) #compile
  thtest2(1,1,1) #force compilation
  println("Testing NO threading")
  #@time
  a, b1, b2, b3 = thtest(1000, 1000, 1000)
  println("testing THREADING")
  #@time
  b, c1, c2, c3 = thtest2(1000, 1000, 1000)
  println("No threading count: ", a)
  println("No threading total: ", a.a + a.b + a.c)
  println("loop counts: ", b1, " ", b2, " ", b3)
  println("Threading count: ", b)
  println("Threading total: ", b.a + b.b + b.c)
  println("Loop counts: ", c1, " ", c2, " ", c3)

end


function test2()
  println("version 2")
  println("compile")
  V1(10)
  V2(10)
  println("Testing NO threads: ")
  @time V1(100000);
  println("Testing THREADS: ")
  @time V2(100000);
end

function test3()
  println("version 3")
  println("compile")
  F1(1)
  F2(1)
  println("testing")
  println("NO threads")
  @time F1(1000);
  println("THREADS")
  @time F2(1000)
end


#doit();
#test2();
#test3();






####



"""
`function NOTHREADChoiceVector(utils::Array{Float64,1},fids::Array{Float64,1},x::Int64)::Array{Float64,1}`
Just for comparison, this one runs WITHOUT threading and takes
about 3-4 times as long for 100_000 entries.  That doesn't mean that
this is the best way you could do this task in a no-threading way.
"""
function NOTHREADChoiceVector(utils::Array{Float64,1},
              fids::Array{Float64,1},
              x::Int64)::Array{Float64,1}
  outp::Array{Float64,1} = zeros(Float64, x)
  temparry::Array{Float64, 1} = zeros(fids)
  for i = 1:x
    outp[i] = UMap(utils, fids, temparry)
  end
  return outp::Array{Float64,1}
end

"""
I think the reason threading doesn't do anything here is that there's only one copy of the dict in memory, so
you can't actually divide up that part of the process across threads.
"""
function FAKETHREADFIDCounter(chosen::Array{Float64,1}, fids::Array{Float64,1})
  outp::Dict{Float64, Int64} = Dict{Float64, Int64}()
  for el in fids
    outp[el] = 0
  end
  Threads.@threads for c in 1:size(chosen,1)
    # this probably doesn't work because there is only one copy of the dictionary.  How can it be mapped
    # across threads?
    outp[chosen[c]] += 1
  end
  return outp
end
