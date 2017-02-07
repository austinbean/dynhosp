

type threadtest
  a::Int64
  b::Int64
  c::Int64
end

#=
After importing the project module.
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);
GenPChoices(patients, Texas);
NB: DicttoVec(patients.zips[xxxx].pdetutils)

=#

function PCh(p::patientcount, )
#TODO - THIS is the part I want to get right.  This can be mapped.
# can the allocation be avoided somehow?
  for k = 1:pats.zips[zipcode].mpatients.count385
    outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count385 += 1
  end


end






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

"""
`function UMap(utils::Array{Float64,2}, fids::Array{Int64,2}, temparr::Array{Float64,2})`

This function takes:
- Array of Utils
- Array of Fids
- Temporary Array

Computes utility + random component, maps out corresponding FID.

testing:
ut = [0.1, 0.2, 0.3, 0.4, 0.5];
fi = [111.0, 222.0, 333.0, 444.0, 19.0];
ta = [0.0, 0.0, 0.0, 0.0, 0.0];
UMap(ut, fi, ta)
"""
function UMap(utils::Array{Float64,1},
              fids::Array{Float64,1},
              temparr::Array{Float64,1};
              dist_μ = 0,
              dist_σ = 1,
              dist_ξ = 0,
              d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))::Float64
  return fids[indmax(utils+rand!(d, temparr))]
end

"""
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


function ChoiceVector(utils::Array{Float64,1},
              fids::Array{Float64,1},
              x::Int64)::Array{Float64,1}
  outp::Array{Float64,1} = zeros(Float64, x)
  temparry::Array{Float64, 1} = zeros(fids)
  Threads.@threads for i = 1:x
    outp[i] = UMap(utils, fids, temparry)
  end
  return outp::Array{Float64,1}
end

function FIDCounter(chosen::Array{Float64,1}, fids::Array{Float64,1})
  outp::Dict{Float64, Int64} = Dict{Float64, Int64}()
  for el in fids
    outp[el] = 0
  end
  for c in chosen
    outp[c] += 1
  end
  return outp
end

# 0.056291 seconds (166.08 k allocations: 15.814 MB)
FIDCounter( ChoiceVector(ab[2,:], ab[1,:], 100000) , ab[1,:] )
# output here is a dictionary.

# What is this going to do?  Take the counts, generate demand by each one,
# then map demand back to patientcounts.
function ZipDemand(z::zip)
  outp::Dict{Float64, patientcount} = Dict{Float64, patientcount}()
  utils::Array{Float64,2} = DicttoVec(z.pdetutils)
  for el in utils[1,:] #TODO - this has to handle medicaid patients too.
    outp[el] =  patientcount(0,0,0,0,0,0,0)
  end
  for k1 in utils[1,:]
    outp[k1].count385 += FIDCounter(ChoiceVector(utils[2,:], utils[1,:], z.ppatients.count385), utils[1,:])[k1]
  end
  for k1 in utils[1,:]
    outp[k1].count386 += FIDCounter(ChoiceVector(utils[2,:], utils[1,:], z.ppatients.count386), utils[1,:])[k1]
  end
  for k1 in utils[1,:]
    outp[k1].count387 += FIDCounter(ChoiceVector(utils[2,:], utils[1,:], z.ppatients.count387), utils[1,:])[k1]
  end
  for k1 in utils[1,:]
    outp[k1].count388 += FIDCounter(ChoiceVector(utils[2,:], utils[1,:], z.ppatients.count388), utils[1,:])[k1]
  end
  for k1 in utils[1,:]
    outp[k1].count389 += FIDCounter(ChoiceVector(utils[2,:], utils[1,:], z.ppatients.count389), utils[1,:])[k1]
  end
  for k1 in utils[1,:]
    outp[k1].count390 += FIDCounter(ChoiceVector(utils[2,:], utils[1,:], z.ppatients.count390), utils[1,:])[k1]
  end
  for k1 in utils[1,:]
    outp[k1].count391 += FIDCounter(ChoiceVector(utils[2,:], utils[1,:], z.ppatients.count391), utils[1,:])[k1]
  end
  return outp
end

function FinalDemand(p::patientcollection, Tex::EntireState)
  outp::Dict{Float64, patientcount} = Dict{Float64, patientcount}()
  outp[0] = patientcount(0,0,0,0,0,0,0)
  for el in keys(Tex.fipsdirectory)
    outp[el] = patientcount(0,0,0,0,0,0,0)
  end
  for z in keys(p.zips)
    for k1 in keys(ZipDemand(p.zips[z]))
      outp[k1] += ZipDemand(p.zips[z])[k1]
    end
  end
  return outp
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
