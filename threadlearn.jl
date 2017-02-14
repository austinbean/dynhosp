



#=
After importing the project module.
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);
GenPChoices(patients, Texas);
NB: DicttoVec(patients.zips[xxxx].pdetutils)

Updated version:
@time ab = OuterSim(3; T1 = 3);
Current iteration 1
Current iteration 2
Current iteration 3
195.408350 seconds (771.19 M allocations: 27.919 GB, 4.49% gc time)

Prior Version:
@time ab = OuterSim(3; T1 = 3);
Current iteration 1
Current iteration 2
Current iteration 3
185.595211 seconds (1.07 G allocations: 48.070 GB, 7.21% gc time)



=#


# Can a type contain a vector of a particular length?
# Not obviously, but you can construct as l1 = lentype(Array{Float64,1}( LENGTH ))
# so the parameter doesn't need to be part of the type.

type lentype2
  elt::Array{Float64,1}
  leng::Int64

  lentype2(elt, leng) = length(elt)==leng ? error("die") : new(elt, leng)
end










# testing pushing vs. indexing.


function pusher(v::Array{Float64,1})
  for i = 1:5000
    push!(v, randn())
  end
end

function writer(v::Array{Float64,1})
  for i = 1:length(v)
    v[i] = randn()
  end
end

function res(v::Array{Float64,1})
  for i = 1:length(v)
    v[i] = 0.0
  end
end

function testp()
  v1 = Array{Float64,1}()
  v2 = Array{Float64,1}(5000)
  pusher(v1)
  writer(v2)
  res(v2)
  for i = 1:10
    v1 = Array{Float64,1}()
    println("Push")
    @time pusher(v1)
    println("Write")
    @time writer(v2)
    res(v2)
  end
end
testp()


function test2()
  v2 = Array{Float64,1}(5000)
  writer(v2)
  res(v2)
  @time for i = 1:50
    v1 = Array{Float64,1}()
    pusher(v1)
  end

  @time for i = 1:50
    writer(v2)
    res(v2)
  end
end
test2() # results suggest that indexing is always faster.  And has no allocations, basically.



# testing immutables:

immutable coeffs
  a::Float64
  b::Float64
  c::Float64
  d::Float64
end

type mutb
  a::Float64
  b::Float64
  c::Float64
  d::Float64
end

type target
  f::Float64
  g::Float64
  h::Float64
  i::Float64
end

m1 = mutb(1.0, 1.0, 1.0, 9.9)
i1 = coeffs(1.0, 1.0, 1.0, 9.9)
t1 = target(2.0, 2.0, 2.0, 2.0)

function f1(m::mutb, t::target)
    return m.a*t.f + m.b*t.g + m.c*t.h + m.d*t.i
end

function f2(c::coeffs, t::target)
  return c.a*t.f + c.b*t.g + c.c*t.h + c.d*t.i
end

f1(m1)
f2(i1)

for i = 1:10
  println("type")
  @time f1(m1)
  println("immutable")
  @time f2(i1)
end



"""
`function UMap(utils::Array{Float64,1}, fids::Array{Int64,1}, temparr::Array{Float64,1})`

This function takes:
- Array of Utils
- Array of Fids
- Temporary Array

Computes utility + random component, maps out corresponding FID.

testing:
#ut = [0.1, 0.2, 0.3, 0.4, 0.5];
#fi = [111, 222, 333, 444, 19];
#ta = [0.0, 0.0, 0.0, 0.0, 0.0];
#UMap(ut, fi, ta)
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
#Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
#patients = NewPatients(Texas);
#DV(patients.zips[78702].pdetutils)
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
`function SingleChoiceVector(utils::Array{Float64,1}, fids::Array{Float64,1}, x::Int64)`
Allocates an array, then uses threading to allocate calls to `UMap` across different threads.  Output
is a Dict{Float64, Int64} of facilities and patient counts.

#Testing on Choice Data:
#Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
#patients = NewPatients(Texas);
#input = ones(Int64, 1650);
#ChoiceVector(patients.zips[78702].pdetutils, input, 50)
"""
function SingleChoiceVector(pd::Dict{Int64, Float64}, ch::Array{Int64,1}, x::Int64)
  fids::Array{Int64,1}, utils::Array{Float64,1} = DV(pd) # this is the biggest source of allocation now.
  # outp::Array{Int64,1} = zeros(Int64, x)
  temparry::Array{Float64, 1} = zeros(utils)
  dt::Dict{Int64, Int64} = Dict()
  for el in fids
    dt[el] = 0
  end
  UseThreads(ch, fids, utils, temparry, x)
  for i = 1:x
    dt[ch[i]] += 1 # It seems like I can't get a key error, but I'm not 100% sure of that.
  end
  ResVec(ch) #reset the vector.
  return dt::Dict{Int64, Int64}
end

"""
`function NewHospDict(Tex::EntireState)`
Associates with each FID a zero patientcount.
The purpose is to make a dictionary to hold the output of ChoiceVector.
"""
function NewHospDict(Tex::EntireState)
  dt::Dict{Int64, patientcount} = Dict()
  for el in keys(Tex.fipsdirectory)
    dt[el] = patientcount(0,0,0,0,0,0,0)
  end
  dt[0] = patientcount(0,0,0,0,0,0,0)
  return dt
end

"""
`function PatientsClean(dt::Dict{Int64, patientcount})`
resets the values of all patientcounts to 0.  This allows re-use of
the dictionary.  Maybe the allocation could be even lower by manually
setting all fields back to 0.  But whatever.
"""
function PatientsClean(dt::Dict{Int64, patientcount})
  for el in keys(dt)
    dt[el] = patientcount(0,0,0,0,0,0,0)
  end
end


"""
`function ChoiceVector(pd::Dict{Int64, Float64},dt::Dict{Int64, patientcount},ch::Array{Int64,1},x::patientcount)`
This should operate in-place on the dictionary dt.
The dictionary has an entry for every hospital.
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);
dic1 = NewHospDict(Texas);
inpt = ones(Int64, 1550); # largest group is 1511
ChoiceVector2(patients.zips[78759].pdetutils, dic1, inpt, patients.zips[78759].ppatients)
REMEMBER TO TURN ON THREADING.
"""
function ChoiceVector(pd::Dict{Int64, Float64},
                      dt::Dict{Int64, patientcount},
                      ch::Array{Int64,1},
                      x::patientcount)
  fids::Array{Int64,1}, utils::Array{Float64,1} = DV(pd)
  temparry::Array{Float64, 1} = zeros(utils)
  for (loc, num) in enumerate(x)
    UseThreads(ch, fids, utils, temparry, num)
    if loc == 1
      for i = 1:num
        dt[ch[i]].count385 += 1
      end
    elseif loc==2
      for i = 1:num
        dt[ch[i]].count386 += 1
      end
    elseif loc==3
      for i = 1:num
        dt[ch[i]].count387 += 1
      end
    elseif loc==4
      for i = 1:num
        dt[ch[i]].count388 += 1
      end
    elseif loc==5
      for i = 1:num
        dt[ch[i]].count389 += 1
      end
    elseif loc==6
      for i = 1:num
        dt[ch[i]].count390 += 1
      end
    elseif loc==7
      for i = 1:num
        dt[ch[i]].count391 += 1
      end
    end
    ResVec(ch) #reset the vector.
  end
end

"""
`function UseThreads(inpt::Array{Int64,1},fids::Array{Int64,1}, utils::Array{Float64,1}, temparry::Array{Float64, 1})`
The sole purpose of this is to avoid #15276, which generates an ambiguity in the type of the arrays in `ChoiceVector`.
There is no ambiguity in this one.
https://github.com/JuliaLang/julia/issues/15276
Workaround, maybe?
See this:
https://github.com/yuyichao/explore/blob/8d52fb6caa745a658f2c9bbffd3b0f0fe4a2cc48/julia/issue-17395/scale.jl#L21
"""
function UseThreads(inpt::Array{Int64,1},fids::Array{Int64,1},utils::Array{Float64,1},temparry::Array{Float64, 1}, x::Int64)
  Threads.@threads for i = 1:x
    inpt[i] = UMap(utils, fids, temparry)
  end
end


"""
`function ResVec(v::Array{Int64,1})`
This function resets the vector to its original state.  This enables it to be
re-used so it does not require so many allocations.
"""
function ResVec(v::Array{Int64, 1})
  for i = 1:length(v)
    v[i] = 1
  end
end

"""
`function GenPChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})`
Returns a dictionary of private demands for the whole state.
Takes as input a patientcollection, a dict{Int64, patientcount} and a re-usable array v.
"""
function GenPChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})
  for k in keys(p.zips)
    ChoiceVector(p.zips[k].pdetutils, d, v, p.zips[k].ppatients)
  end
end

"""
`function GenMChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})`
Returns a dictionary of private demands for the whole state.
Takes as input a patientcollection, a dict{Int64, patientcount} and a re-usable array v.
"""
function GenMChoices(p::patientcollection, d::Dict{Int64, patientcount}, v::Array{Int64,1})
  for k in keys(p.zips)
    ChoiceVector(p.zips[k].pdetutils, d, v, p.zips[k].ppatients)
  end
end








function testit()
  Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
  patients = NewPatients(Texas);
  dic1 = NewHospDict(Texas)
  inp1 = zeros(Int64, 1550)
  StateDemand(patients, dic1, inp1)
  println("Timing")
  @time StateDemand(patients, dic1, inp1)
  PatientsClean(dic1)
  for i = 1:5
    println("Old")
    @time GenPChoices(patients, Texas)
    println("New")
    @time StateDemand(patients, dic1, inp1)
    PatientsClean(dic1)
  end

end

using Distributions

testit()



"""
`function DictCombine(outp, arg...)`
This function takes an arbitrary number of dictionaries and adds them to outp.
Dicts must be Dict{Int64, Int64}.

#D1 = Dict(1 => 3, 2 => 4, 3 => 16, 4 => 12)
#D2 = Dict( 6=> 36, 8=>64 , 10=>100 , 12=> 144, 14=>1414 )
#D3 = Dict( 7=> 49, 9=> 81, 11=> 121, 13=> 169, 15=> 225)
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

function DoIt()
  #D::Dict{Int64, Float64} = Dict(4536048 => 0.17909, 4536337 => -0.00907679, 4536253 => -0.00466823, 4916029 => 0.155892, 0 => 0.0, 4530190 => 1.45951, 2093151 => 2.12427, 4916068 => -0.0048307, 4536338 => 0.346297, 4530200 => 1.64229, 4530170 => 1.5257)
  ChoiceVector(Dict(4536048 => 0.17909, 4536337 => -0.00907679, 4536253 => -0.00466823, 4916029 => 0.155892, 0 => 0.0, 4530190 => 1.45951, 2093151 => 2.12427, 4916068 => -0.0048307, 4536338 => 0.346297, 4530200 => 1.64229, 4530170 => 1.5257), 1)
  ChoiceVector(Dict(4536048 => 0.17909, 4536337 => -0.00907679, 4536253 => -0.00466823, 4916029 => 0.155892, 0 => 0.0, 4530190 => 1.45951, 2093151 => 2.12427, 4916068 => -0.0048307, 4536338 => 0.346297, 4530200 => 1.64229, 4530170 => 1.5257), 100_000);
end

#DoIt()


function UMapAgg(utils::Array{Float64,1},
              fids::Array{Int64,1},
              temparr::Array{Float64,1},
              x::Int64) # NB - could change this int to be an array or a patientcount, etc.
  out1::Dict{Int64,Int64} = Dict()
  for el in fids
    out1[el] = 0
  end
  for i = 1:x
    out1[UMap(utils, fids, temparr)] += 1 #recall UMap returns a fid.
  end
  return out1
end



# Can threads be called to compute UMap w/ out allocating the vector?
function ThreadDict(pd::Dict{Int64,Float64}, x::Int64)
  fids::Array{Int64,1}, utils::Array{Float64,1} = DV(pd)
  temparry::Array{Float64, 1} = zeros(fids)
  d1::Dict{Int64,Int64} = Dict()
  for el in keys(pd)
    d1[el] = 0
  end
  # TODO - below doesn't work.   Why not?  Well,
  Threads.@threads for i = 1:x
    println("Thread ", Threads.threadid(), "  ", DictCombine(d1, UMapAgg(utils, fids, temparry, x)))
  end
  #return d1
end

#another idea - make a version of UMap that takes an integer argument and returns a dict.





#=
# Each zip returns a dictionary of {fid, demand} for ONE DRG.
# work with that simple case and then make it worse.
function attempt2(pz::patientcollection)
  outp::Dict{Float64,Int64} = Dict{Float64, Int64}()
  for ky in keys(pz.zips)
    DictCombine(outp, ChoiceVector(pz.zips[ky].pdetutils, pz.zips[ky].ppatients.count391))
  end
  return outp
end
=#


#=
#OLD DEMAND FUNCTIONS

Speed up is on the order of:
Old
  0.307162 seconds (1.08 M allocations: 64.515 MB)
New
  0.175351 seconds (351.06 k allocations: 24.492 MB)
Old
  0.292307 seconds (1.08 M allocations: 64.512 MB)
New
  0.185039 seconds (351.68 k allocations: 24.514 MB)


"""
`GenPChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))`
The patient choice is max \bar{U} + ϵ, but we have \bar{U} from Compute Det Util and we know how many patients there are in
the privately insured category from FillPPatients.  This returns a dict of fids and patient counts, where patient counts are
generated by repeatedly finding max i = 1, ..., N \bar{U}_i + ϵ_i.  Note the corresponding GenMChoices below.

# Testing - to generate the patientcollection
Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);
GenPChoices(patients, Texas);
"""
function GenPChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  outp::Dict{Int64, patientcount} = Dict( j=> patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory) )
  outp[0] = patientcount(0,0,0,0,0,0,0) # create an outside option
  for zipcode in keys(pats.zips)
    if pats.zips[zipcode].pdetutils.count > 0
      utils = hcat([ [k1,pats.zips[zipcode].pdetutils[k1]] for k1 in keys(pats.zips[zipcode].pdetutils)]...)
      temparr = zeros(size(utils, 2))
      for k = 1:pats.zips[zipcode].ppatients.count385
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count385 += 1
      end
      for k = 1:pats.zips[zipcode].ppatients.count386
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count386 += 1
      end
      for k = 1:pats.zips[zipcode].ppatients.count387
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count387 += 1
      end
      for i=1:pats.zips[zipcode].ppatients.count388
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count388 += 1
      end
      for i = 1:pats.zips[zipcode].ppatients.count389
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count389 += 1
      end
      for i=1:pats.zips[zipcode].ppatients.count390
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count390 += 1
      end
      for i = 1:pats.zips[zipcode].ppatients.count391
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count391 += 1
      end
    end
  end
  return outp
end




"""
`GenMChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))`
The patient choice is max \bar{U} + ϵ, but we have \bar{U} from Compute Det Util and we know how many patients there are in
the Medicaid category from FillMPatients.  This returns a dict of fids and patient counts, where patient counts are
generated by repeatedly finding max i = 1, ..., N \bar{U}_i + ϵ_i.  Note the corresponding GenPChoices above.
This function does permit patients to choose the outside option.

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
patients = NewPatients(Texas);
GenMChoices(patients, Texas)
"""
function GenMChoices(pats::patientcollection, Tex::EntireState; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  # Is the next line really slow?
   outp::Dict{Int64, patientcount} = Dict( j=> patientcount(0, 0, 0, 0, 0, 0, 0) for j in keys(Tex.fipsdirectory) )
  #outp::Dict{Int64, patientcount} = PatientCountOut(Tex)
  outp[0] = patientcount(0, 0, 0, 0, 0, 0, 0) # adding a zero entry - patients are permitted to choose the outside option.
  for zipcode in keys(pats.zips)
    if pats.zips[zipcode].mdetutils.count > 0
      utils = DicttoVec(pats.zips[zipcode].mdetutils) # this creates a 2X(num hospitals) matrix - top row is fids, bottom is deterministic utilities.
      temparr = zeros(size(utils, 2))
      for k = 1:pats.zips[zipcode].mpatients.count385
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count385 += 1
      end
      for k = 1:pats.zips[zipcode].mpatients.count386
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count386 += 1
      end
      for k = 1:pats.zips[zipcode].mpatients.count387
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count387 += 1
      end
      for i=1:pats.zips[zipcode].mpatients.count388
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count388 += 1
      end
      for i = 1:pats.zips[zipcode].mpatients.count389
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count389 += 1
      end
      for i=1:pats.zips[zipcode].mpatients.count390
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count390 += 1
      end
      for i = 1:pats.zips[zipcode].mpatients.count391
        outp[utils[1,indmax(utils[2,:] + rand!(d, temparr))]].count391 += 1
      end
    end
  end
  return outp
end




=#
