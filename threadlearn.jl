



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
`function ChoiceVector(utils::Array{Float64,1}, fids::Array{Float64,1}, x::Int64)`
Allocates an array, then uses threading to allocate calls to `UMap` across different threads.  Output
is a Dict{Float64, Int64} of facilities and patient counts.

#Testing on Choice Data:
#Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists);
#patients = NewPatients(Texas);

#ChoiceVector(patients.zips[78702].pdetutils, 50)
#NB - there is a performance issue in that the type of many of these is read as Core.Box, but this is #15276 or something.
"""
function ChoiceVector(pd::Dict{Int64, Float64}, x::Int64)
  fids::Array{Int64,1}, utils::Array{Float64,1} = DV(pd)
  outp::Array{Int64,1} = zeros(Int64, x)
  temparry::Array{Float64, 1} = zeros(fids)
  dt::Dict{Int64, Int64} = Dict()
  for el in fids
    dt[el] = 0
  end
  # Uncomment when 15276 is fixed.
  # Threads.@threads for i = 1:x
  #   outp[i] = UMap(utils, fids, temparry)
  # end
  UseThreads(outp, fids, utils, temparry)
  for i = 1:x
    dt[outp[i]] += 1
  end
  return dt::Dict{Int64, Int64}
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
function UseThreads(inpt::Array{Int64,1},fids::Array{Int64,1}, utils::Array{Float64,1}, temparry::Array{Float64, 1})
  Threads.@threads for i = 1:length(inpt)
    inpt[i] = UMap(utils, fids, temparry)
  end
end



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

using Distributions
DoIt()





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
