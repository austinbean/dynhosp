

type threadtest
  a::Int64
  b::Int64
  c::Int64
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

function mapper(utils::Array{Float64,2}, temparr::Array{Float64,2})
  return indmax(utils+randn!(temparr))
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
  return a
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
  return a
end




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
  V1(1)
  V2(1)
  println("Testing NO threads: ")
  @time V1(1000);
  println("Testing THREADS: ")
  @time V2(1000);
end



#doit()
test2()
