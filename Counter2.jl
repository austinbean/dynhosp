# Now solve the dynamic model itself.


# Can put these into a market and then into a state.

# use the type neighbors !

using Distributions
using DataFrames

"""
`type nlrec`
- aw::Dict{Int64,Float64}
- psi::Array{Float64,2}
- counter::Dict{Int64, Int64}
First element stores {Action, Continuation Value Approximations}.  Second stores probabilities.  Third is an {Action, Hits} counter.
Stored in a dictionary under a (neighbors, level) tuple-type key.
"""
type nlrec
  aw::Dict{Int64,Float64}
  psi::Array{Float64,2}
  counter::Dict{Int64, Int64}
end

"""
`type hitcount`
- conf::neighbors
- visits::Dict{Int64, Int64}
"""
type hitcount # this type will record visits to a state-action pair
  conf::neighbors
  visits::Dict{Int64,Int64}
end


"""
`type history`
-  path::Dict{neighbors, hitcount}
-  totalcount::Int64
"""
type history
  path::Dict{neighbors, hitcount}
  totalcount::Int64 # records total number of iterations
end


"""
`type shortrec<:ProjectModule.Fac`
-  fid::Int64
-  lat::Float64
-  long::Float64
-  level::Int64
-  truelevel::Int64
-  beds::Int64
-  ns::neighbors
-  choices::Array{Int64, 2}
-  chprobs::WeightVec
-  tbu::Bool
"""
type shortrec<:ProjectModule.Fac
  # needs to take location.  must measure distances.  need to update all of these every time, probably.
  fid::Int64
  lat::Float64
  long::Float64
  level::Int64
  truelevel::Int64
  beds::Int64
  ns::neighbors
  choices::Array{Int64, 2}
  chprobs::WeightVec
  tbu::Bool
end

#NB:  this array needs to include the WTP for each facility too!
"""
`type cpats`
-  zp::Int64
-  lat::Float64
-  long::Float64
-  putils::Array{Float64,2}
-  mutils::Array{Float64,2}
-  pwtp::Array{Float64,2}
-  facs::Array{shortrec,1}
-  pcounts::patientcount
-  mcounts::patientcount
"""
type cpats
  zp::Int64
  lat::Float64
  long::Float64
  putils::Array{Float64,2}
  mutils::Array{Float64,2}
  pwtp::Array{Float64,2}
  facs::Array{shortrec,1}
  pcounts::patientcount
  mcounts::patientcount
end


"""
`type cmkt`
-  fid::Int64
-  m::Array{cpats,1}
"""
type cmkt
  fid::Int64
  m::Array{cpats,1}
end


#NB: When the firm exits, can probably restart from the beginning, but keeping the elements in "visited".  We can keep approximating them.
"""
`type simh<:ProjectModule.Fac`
-  fid::Int64
-  lat::Float64
-  long::Float64
-  level::Int64
-  actual::Int64
-  beds::Int64
-  cns::neighbors # must know what current neighbors look like.
-  visited::Dict{nl, nlrec} #possible given "isequal" and "hash" extended for "neighbors"
-  ns::Array{shortrec, 1}
-  mk::cmkt # putting the cmkt into the simh record itself.
-  exit::Bool
-  tbu::Bool
"""
type simh<:ProjectModule.Fac
  fid::Int64
  lat::Float64
  long::Float64
  level::Int64
  actual::Int64
  beds::Int64
  cns::neighbors # must know what current neighbors look like.
  visited::Dict{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, nlrec} #change key to tuple of int64 composed of neighbors and level.
  ns::Array{shortrec, 1}
  mk::cmkt # putting the cmkt into the simh record itself.
  exit::Bool
  tbu::Bool
end

"""
`type DynState` # this should hold a collection of ALL of the hospitals, for the approximation.
-  all::Array{simh, 1}
"""
type DynState # this should hold a collection of ALL of the hospitals, for the approximation.
  all::Array{simh, 1}
end


"""
`KeyCreate(n::neighbors, l::Int64)`
takes neighbors and level and returns a tuple - this is a better way (hopefully) to make the keys for the nlrec.
"""
function KeyCreate(n::neighbors, l::Int64)
  return (n.level105, n.level205, n.level305, n.level1515, n.level2515, n.level3515, n.level11525, n.level21525, n.level31525, l)
end



"""
`DynStateCreate(Tex::EntireState, Tex2::EntireState, p::patientcollection )`
Create the dynamic records from the existing state, don't bother doing it from scratch.
And these don't have to be organized by zip.
Use the EntireState from the equilibrium simulation, not the first counterfactual.
Make using
TexasEq = MakeNew(ProjectModule.fips, ProjectModule.data05);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.data05);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);

Note that the function takes TWO EntireState arguments.  This is super dumb, but
only one of them (containing hospital types) has the bed counts.
"""
#TODO - need to make sure an outside option is added to each zipcode record.
function DynStateCreate( Tex::EntireState, Tex2::EntireState, p::patientcollection )
  outp = DynState(Array{simh,1}())
  for k1 in keys(Tex.mkts)
    for hk in keys(Tex.mkts[k1].collection)
      newsimh = simh(Tex.mkts[k1].collection[hk].fid,
                     Tex.mkts[k1].collection[hk].lat,
                     Tex.mkts[k1].collection[hk].long,
                     Tex.mkts[k1].collection[hk].level,
                     Tex.mkts[k1].collection[hk].level,
                     convert(Int64, Tex2.mkts[k1].collection[hk].bedcount),
                     neighbors(0,0,0,0,0,0,0,0,0),
                     Dict{Tuple{Int64}, nlrec}(),
                     Array{shortrec,1}(),
                     DynPatients(p, Tex.mkts[k1].collection[hk].fid), # should create the patient collection as a subelement of the hospital record.
                     false,
                     false)
      for k2 in keys(Tex.mkts)
        for hk2 in keys(Tex.mkts[k2].collection)
          if hk2 != hk
            d1 = distance(Tex.mkts[k1].collection[hk].lat, Tex.mkts[k1].collection[hk].long , Tex.mkts[k2].collection[hk2].lat , Tex.mkts[k2].collection[hk2].long)
            if d1 < 25 #check distance
              push!(newsimh.ns, shortrec(Tex.mkts[k2].collection[hk2].fid,
                                         Tex.mkts[k2].collection[hk2].lat,
                                         Tex.mkts[k2].collection[hk2].long,
                                         Tex.mkts[k2].collection[hk2].level,
                                         Tex.mkts[k2].collection[hk2].level,
                                         convert(Int64, Tex2.mkts[k2].collection[hk2].bedcount),
                                         neighbors(0,0,0,0,0,0,0,0,0),
                                         ChoicesAvailable(Tex.mkts[k2].collection[hk2]),
                                         Tex.mkts[k2].collection[hk2].chprobability,
                                         false))
              if (d1>0)&(d1<5)
                if Tex.mkts[k2].collection[hk2].level == 1
                  newsimh.cns.level105 += 1
                elseif Tex.mkts[k2].collection[hk2].level == 2
                  newsimh.cns.level205 += 1
                else # equals 3
                  newsimh.cns.level305 += 1
                end
              elseif (d1>=5)&(d1<15)
                if Tex.mkts[k2].collection[hk2].level == 1
                  newsimh.cns.level1515 += 1
                elseif Tex.mkts[k2].collection[hk2].level == 2
                  newsimh.cns.level2515 += 1
                else # equals 3
                  newsimh.cns.level3515 += 1
                end
              else # d1 >= 15 & d1 < 25
                if Tex.mkts[k2].collection[hk2].level == 1
                  newsimh.cns.level11525 += 1
                elseif Tex.mkts[k2].collection[hk2].level == 2
                  newsimh.cns.level21525 += 1
                else # equals 3
                  newsimh.cns.level31525 += 1
                end
              end
            end
          end
        end
      end
      # This functionality now duplicated by FixNN(s::simh) below.
      for el in newsimh.ns
        for el2 in newsimh.ns
          if el.fid != el2.fid
            if (distance(el.lat, el.long, el2.lat, el2.long)>0)&(distance(el.lat, el.long, el2.lat, el2.long)<5)
              if el2.level == 1
                el.ns.level105+=1
              elseif el2.level == 2
                  el.ns.level205+=1
              else # equals 3
                  el.ns.level305+=1
              end
            elseif (distance(el.lat, el.long, el2.lat, el2.long)>=5)&(distance(el.lat, el.long, el2.lat, el2.long)<15)
              if el2.level == 1
                  el.ns.level1515+=1
              elseif el2.level == 2
                  el.ns.level2515+=1
              else # equals 3
                  el.ns.level3515+=1
              end
            elseif (distance(el.lat, el.long, el2.lat, el2.long)>=15)&(distance(el.lat, el.long, el2.lat, el2.long)<=25)
              if el2.level == 1
                  el.ns.level11525+=1
              elseif el2.level == 2
                  el.ns.level21525+=1
              else # equals 3
                  el.ns.level31525+=1
              end
            end
          end
        end
      end
      for el in newsimh.ns       # these are shortrecs
        for zp in newsimh.mk.m   # these are cpats
          if in(el.fid, zp.putils[1,:])
            push!(zp.facs, el)   # This should push all of the shortrec types
          end
        end
      end
      push!(outp.all, newsimh)
    end
  end
  return outp
end


"""
`DynPatients(p::patientcollection, f::fid)`
This will take a collection of patients and create a cmkt, which is a vector of
cpats.  That is, it will take every zip for which `f` is an option, then create
the collection of patients for those zips.  Note that `f` is a FID for a hospital.
"""
function DynPatients(p::patientcollection, f::Int64 )
  outp::cmkt = cmkt(f, Array{cpats,1}())
  zpc = PatientFind(p, f) # finds the zip codes
  for el in zpc # TODO: add the outside option here.
    push!(outp.m, cpats(el,
                  p.zips[el].lat,
                  p.zips[el].long,
                  DetUtils(p.zips[el]; switch = false),
                  DetUtils(p.zips[el]; switch = true),
                  vcat(transpose(DetUtils(p.zips[el]; switch = false)[1,:]), transpose(CounterWTP(DetUtils(p.zips[el]; switch = false)[2,:]))), #NB: bottom row only.
                  Array{shortrec,1}(),
                  p.zips[el].ppatients,
                  p.zips[el].mpatients ) ) #note - this is *not* a copy
  end
  return outp
end


"""
`DetUtils(z::zip)`
Returns a 2 x N array of the (fid, utility) pairs from the zipcode z
Top row will be the collection of fids
and the bottom row will be the utilities themselves.
Outside option is added here as [0, 0].
"""
function DetUtils(z::zip; switch::Bool = false)
  if switch
    return  hcat(hcat([[k, z.pdetutils[k]] for k in keys(z.pdetutils)]...),[0,0])
  else
    return hcat(hcat([[k, z.mdetutils[k]] for k in keys(z.mdetutils)]...),[0,0])
  end
end



"""
`CounterWTP(ar::Array{Float64,2})`
Will take the output of `DetUtils(z::zipc; switch)` and compute the WTP.
This only needs to be done for the private patients!
"""
function CounterWTP(ar::Array{Float64})
  denom::Float64 = 0.0
  for i =1:maximum(size(ar)) #might be a problem for one choice.
    # try exp(ar[i])
    # catch y
    #   println(y)
    #   if y == DomainError
    #     println(ar[i])
    #   end
    # end
    denom += (ar[i] = exp(ar[i]))
  end
  return log(map(x->(1/(1-x)), ar./denom))
end

"""
`DSim(c::cmkt, f::Int64)`
This is going to take the collection of patients and the fid and figure out
how many people choose it.  It's ok fast, but not really fast.
"""
function DSim(c::cmkt, f::Int64; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  pcount::patientcount = patientcount(0,0,0,0,0,0,0)
  mcount::patientcount = patientcount(0,0,0,0,0,0,0)
  for el in c.m
    # NB: here is a parallel opportunity, maybe?  Sum across the zip codes across cores?  Or a threading opportunity?
    # Threading possible *per loop* perhaps possible using the Threads.@threads annotation on each one?
    siz1 = size(el.putils[2,:],1) #siz1 and siz2 should always be the same.
    siz2 = size(el.mutils[2,:],1)
    for i = 1:el.pcounts.count385
      if el.putils[1,indmax( el.putils[2,:] +rand(d, siz1) )] == f
        pcount.count385 += 1
      end
    end
    for i = 1:el.pcounts.count386
      if el.putils[1,indmax( el.putils[2,:] +rand(d, siz1) )] == f
        pcount.count386 += 1
      end
    end
    for i = 1:el.pcounts.count387
      if el.putils[1,indmax( el.putils[2,:] +rand(d, siz1) )] == f
        pcount.count387 += 1
      end
    end
    for i = 1:el.pcounts.count388
      if el.putils[1,indmax( el.putils[2,:] +rand(d, siz1) )] == f
        pcount.count388 += 1
      end
    end
    for i = 1:el.pcounts.count389
      if el.putils[1,indmax( el.putils[2,:] +rand(d, siz1) )] == f
        pcount.count389 += 1
      end
    end
    for i = 1:el.pcounts.count390
      if el.putils[1,indmax( el.putils[2,:] +rand(d, siz1) )] == f
        pcount.count390 += 1
      end
    end
    for i = 1:el.pcounts.count391
      if el.putils[1,indmax( el.putils[2,:] +rand(d, siz1) )] == f
        pcount.count391 += 1
      end
    end
    #NB: Medicaid patients here:
    for i = 1:el.mcounts.count385
      if el.mutils[1,indmax( el.mutils[2,:] +rand(d, siz2) )] == f
        mcount.count385 += 1
      end
    end
    for i = 1:el.mcounts.count386
      if el.mutils[1,indmax( el.mutils[2,:] +rand(d, siz2) )] == f
        mcount.count386 += 1
      end
    end
    for i = 1:el.mcounts.count387
      if el.mutils[1,indmax( el.mutils[2,:] +rand(d, siz2) )] == f
        mcount.count387 += 1
      end
    end
    for i = 1:el.mcounts.count388
      if el.mutils[1,indmax( el.mutils[2,:] +rand(d, siz2) )] == f
        mcount.count388 += 1
      end
    end
    for i = 1:el.mcounts.count389
      if el.mutils[1,indmax( el.mutils[2,:] +rand(d, siz2) )] == f
        mcount.count389 += 1
      end
    end
    for i = 1:el.mcounts.count390
      if el.mutils[1,indmax( el.mutils[2,:] +rand(d, siz2) )] == f
        mcount.count390 += 1
      end
    end
    for i = 1:el.mcounts.count391
      if el.mutils[1,indmax( el.mutils[2,:] +rand(d, siz2) )] == f
        mcount.count391 += 1
      end
    end
  end
  return pcount, mcount
end



"""
`UpdateDUtil(h::simh)`
When something in the market changes, the utility must be updated in all zip codes
for which it can be chosen.  This calls the function HUtil to do the actual update.

"""
function UpdateDUtil(h::simh)
  for el in h.mk.m
    for sr in el.facs
      if sr.tbu # if this is true, the hospital needs updating
        el.putils[2, findin(el.putils[1,:], sr.fid)] = HUtil(el, sr, true)
        el.mutils[2, findin(el.mutils[1,:], sr.fid)] = HUtil(el, sr, false)
      end
      sr.tbu = false # reset.
    end
    if h.tbu # if this is true, the main hosp needs updating
      el.putils[2, findin(el.putils[1,:], h.fid)] = HUtil(el, h, true)
      el.mutils[2, findin(el.mutils[1,:], h.fid)] = HUtil(el, h, false)
    end
  end
  h.tbu = false; # reset the main hospital.
end



"""
`HUtil{T<:Fac}( c::cpats, sr::T, p_or_m::Bool; pparameters::ProjectModule.coefficients = ..., mparameters::ProjectModule.coefficients  = ...)`
Computes the deterministic component of the utility for a hospital - the goal is for this to update the entire cpat list of facilities.
And I want to do both - I don't want to update either private or medicaid.
"""

function HUtil{T<:ProjectModule.Fac}(c::cpats, sr::T, p_or_m::Bool;
               mcoeffs::ProjectModule.coefficients =  coefficients(ProjectModule.medicaiddistance_c, ProjectModule.medicaiddistsq_c, ProjectModule.medicaidneoint_c, ProjectModule.medicaidsoloint_c, ProjectModule.medicaiddistbed_c, ProjectModule.medicaidclosest_c),
               pcoeffs::ProjectModule.coefficients = coefficients(ProjectModule.privatedistance_c, ProjectModule.privatedistsq_c, ProjectModule.privateneoint_c, ProjectModule.privatesoloint_c, ProjectModule.privatedistbed_c, ProjectModule.privateclosest_c),)
  d::Float64 = distance(c.lat, c.long, sr.lat, sr.long)
  if p_or_m # if TRUE private
    if sr.level == 1
      return pcoeffs.distance*d+pcoeffs.distsq*(d^2)+pcoeffs.distbed*(sr.beds*d/100)+pcoeffs.closest*(0)
    elseif sr.level == 2
      return pcoeffs.distance*d+pcoeffs.distsq*(d^2)+pcoeffs.distbed*(sr.beds*d/100)+pcoeffs.closest*(0)+pcoeffs.inter
    else #sr.level == 3
      return pcoeffs.distance*d+pcoeffs.distsq*(d^2)+pcoeffs.distbed*(sr.beds*d/100)+pcoeffs.closest*(0)+pcoeffs.inten
    end
  else
    if sr.level == 1
      return mcoeffs.distance*d+mcoeffs.distsq*(d^2)+mcoeffs.distbed*(sr.beds*d/100)+mcoeffs.closest*(0)
    elseif sr.level == 2
      return mcoeffs.distance*d+mcoeffs.distsq*(d^2)+mcoeffs.distbed*(sr.beds*d/100)+mcoeffs.closest*(0)+mcoeffs.inter
    else #sr.level == 3
      return mcoeffs.distance*d+mcoeffs.distsq*(d^2)+mcoeffs.distbed*(sr.beds*d/100)+mcoeffs.closest*(0)+mcoeffs.inten
    end
  end
end



"""
`UpdateCheck(h::simh)`
Takes a simh record `h` and checks to see whether anyfacility needs updating.
Note that doing it this way seems to work and is much faster than the function
`reduce(&, [el.tbu for el in h.ns])` which allocates the array.
"""
function UpdateCheck(h::simh)
  b::Bool = h.tbu
  if !b
    for el in h.ns
      if el.tbu
        b = true
        break
      end
    end
  end
  return b
end





"""
`NCheck(d::DynState, e::EntireState)`
This function just prints out the records of hospitals which have one or two neighbors to check them against
the same record for the same hospital in the EntireState record type.    Solely for debugging purposes.
"""
function NCheck(d::DynState, e::EntireState)
  for el in d.all
    if size(el.ns,1)<2
      println("*************************")
      println(el.fid)
      println(el.ns)
      println(e.mkts[e.fipsdirectory[el.fid]].collection[el.fid].hood)
    end
  end
end


"""
`GetProb(s::simh)`
Take the simh hospital record, get some action choice by the other firms, update their choices and records.
Return the state, probably.
NB: this does not include the possibility of entry so far.
"""
function GetProb(s::simh)
  #TODO - eventually this should allow entry too.
  for el in s.ns
    if el.level != -999
      action = sample(ChoicesAvailable(el), el.chprobs )                            # Take the action
      di = distance(el.lat, el.long, s.lat, s.long)
      if action == 1
        el.level = 3
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        # NB: this is not working exactly correctly since some of these levels are negative!
        el.chprobs = WeightVec(vec(logitest((0,1), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        if (di>0)&(di<5)
          s.cns.level305 += 1
          s.cns.level105 -= 1
        elseif (di>=5)&(di<15)
          s.cns.level3515 += 1
          s.cns.level1515 -= 1
        elseif (di>=15)&(di<25)
          s.cns.level31525 += 1
          s.cns.level11525 -= 1
        else #nothing here
          #nothing
        end
        el.tbu = true
      elseif action == 2
        el.level = 2
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((1,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        if (di>0)&(di<5)
          s.cns.level205 += 1
          s.cns.level105 -= 1
        elseif (di>=5)&(di<15)
          s.cns.level2515 += 1
          s.cns.level1515 -= 1
        elseif (di>=15)&(di<25)
          s.cns.level21525 += 1
          s.cns.level11525 -= 1
        else #nothing here
          #nothing
        end
        el.tbu = true
      elseif action == 3
        el.level = 2
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((1,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        if (di>0)&(di<5)
          s.cns.level205 += 1
          s.cns.level305 -= 1
        elseif (di>=5)&(di<15)
          s.cns.level2515 += 1
          s.cns.level3515 -= 1
        elseif (di>=15)&(di<25)
          s.cns.level21525 += 1
          s.cns.level31525 -= 1
        else #nothing here
          #nothing
        end
        el.tbu = true
      elseif action == 4
        el.level = 1
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((0,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        if (di>0)&(di<5)
          s.cns.level105 += 1
          s.cns.level305 -= 1
        elseif (di>=5)&(di<15)
          s.cns.level1515 += 1
          s.cns.level3515 -= 1
        elseif (di>=15)&(di<25)
          s.cns.level11525 += 1
          s.cns.level31525 -= 1
        else #nothing here
          #nothing
        end
        el.tbu = true
      elseif action == 5
        el.level = 1
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((0,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        if (di>0)&(di<5)
          s.cns.level105 += 1
          s.cns.level205 -= 1
        elseif (di>=5)&(di<15)
          s.cns.level1515 += 1
          s.cns.level2515 -= 1
        elseif (di>=15)&(di<25)
          s.cns.level11525 += 1
          s.cns.level21525 -= 1
        else #nothing here
          #nothing
        end
        el.tbu = true
      elseif action == 6
        el.level = 3
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((0,1), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        if (di>0)&(di<5)
          s.cns.level305 += 1
          s.cns.level205 -= 1
        elseif (di>=5)&(di<15)
          s.cns.level3515 += 1
          s.cns.level2515 -= 1
        elseif (di>=15)&(di<25)
          s.cns.level31525 += 1
          s.cns.level21525 -= 1
        else #nothing here
          #nothing
        end
        el.tbu = true
      elseif action == 11
        if el.level == 1
          if (di>0)&(di<5)
            s.cns.level105 -= 1
          elseif (di>=5)&&(di<15)
            s.cns.level1515 -= 1
          elseif (di>=15)&&(di<25)
            s.cns.level11525 -= 1
          else #nothing
            #nothing
          end
        elseif el.level == 2
          if (di>0)&(di<5)
            s.cns.level205 -= 1
          elseif (di>=5)&&(di<15)
            s.cns.level2515 -= 1
          elseif (di>=15)&&(di<25)
            s.cns.level21525 -= 1
          else #nothing
            #nothing
          end
        else # level is 3
          if (di>0)&(di<5)
            s.cns.level305 -= 1
          elseif (di>=5)&&(di<15)
            s.cns.level3515 -= 1
          elseif (di>=15)&&(di<25)
            s.cns.level31525 -= 1
          else #nothing
            #nothing
          end
        end
        el.level = -999
        el.tbu = true
      else # action == 10 - "do nothing" (recall that 7,8,9 were entry actions)
        #nothing.
      end
    end
    FixNN(s)   # this will correct the "neighbors of the neighbors"
  end
end



"""
`FixNN(s::simh)`
For the neighbors in the array for s, make sure their `neighbors` type is correctly
counted.
"""
function FixNN(s::simh)
  for el1 in s.ns
    el1.ns = neighbors(0,0,0,0,0,0,0,0,0)
    for el2 in s.ns
      if el1.fid != el2.fid
        dd = distance(el1.lat, el1.long, el2.lat, el2.long)
        if dd < 25 # should always be true
          if  (dd>0)&&(dd<5)
            if el2.level == 1
              el1.ns.level105 += 1
            elseif el2.level == 2
              el1.ns.level205 += 1
            elseif el2.level == 3
              el1.ns.level305 += 1
            else #level == -999 (exited)
              #nothing.
            end
          elseif (dd>=5)&&(dd<15)
            if el2.level == 1
              el1.ns.level1515 += 1
            elseif el2.level == 2
              el1.ns.level2515 += 1
            elseif el2.level == 3
              el1.ns.level3515 += 1
            else #level == -999 (exited)
              #nothing.
            end
          else (dd>=15)&&(dd<25)
            if el2.level == 1
              el1.ns.level11525 += 1
            elseif el2.level == 2
              el1.ns.level21525 += 1
            elseif el2.level == 3
              el1.ns.level31525 += 1
            else #level == -999 (exited)
              #nothing.
            end
          end
        end
      end
    end
  end
end



"""
`PatientFind(s::patientcollection, f::fid)`
Take an entire state type and figure out which zips have that as an option.
Returns a vector of the zips in which it is possible to choose that one hospital
denoted by f.
"""
function PatientFind(s::patientcollection, f::Int64)
  outp = Array{Int64,1}()
  for zipc in keys(s.zips)
    for fid in keys(s.zips[zipc].facilities)
      if (fid == f)&(!in(zipc,outp))
        push!(outp, zipc)
      end
    end
  end
  return outp
end









"""
`GetProbCheck(d::DynState)`
Just simulates the function `GetProb` and tries to figure out whether it's working.
I'd like to know whether the ns actually changes and whether the neighbors get updated...
How am I going to do that?
"""
function GetProbCheck(d::DynState; nsim = 50, fi = ProjectModule.fips, da = ProjectModule.data05)
  d2 = DynStateCreate(MakeNew(fi, da))
  for i = 1:nsim
    for el in d.all
      GetProb(el)
    end
  end
  for el2 in d2.all
    for el in d.all
      if el.fid == el2.fid
        println("********")
        println(typeof(el))
        println(typeof(el2))
        println(el.cns)
        println(el2.cns)
      end
    end
  end
end






"""
`FindWTP(h::simh)`
Get the WTP across all of the zips in the market.  This should search through the
mkt and add up the values
"""
function FindWTP(h::simh)
  fid::Int64 = h.fid
  WTP::Float64 = 0.0
  for el in h.mk.m
    for f in 1:size(el.pwtp,2)
      if el.pwtp[1,f] == fid
        WTP += el.pwtp[2,f]
      end
    end
  end
  return WTP
end

"""
`CatchWTP(h::simh)`
When WTP is infinite, what went wrong?
"""
function CatchWTP(h::simh; print::Bool = false)
  outp::Array{Int64,1} = Array{Int64,1}()
  for el in h.mk.m
    for f in 1:size(el.pwtp, 2)
      if el.pwtp[2,f] == Inf
        push!(outp, h.fid)
        if print
          println(h.fid)
          println(el)
          println(el.pwtp)
        end
      end
    end
  end
  return outp
end

"""
`CatchWTPAll(dyn::DynState)`
get the fids of all of the firms with infinite wtp.
"""
function CatchWTPAll(dyn::DynState)
  outp::Array{Int64,1} = Array{Int64,1}()
  outp2::Array{Int64,1} = Array{Int64,1}()
  for el in 1:maximum(size(dyn.all))
    for el2 in CatchWTP(dyn.all[el])
      push!(outp, el2)
      push!(outp2, el)
    end
  end
  return outp, unique(outp2)
end



"""
`SinglePay(s::simh, mpats::ProjectModule.patientcount, ppats::ProjectModule.patientcount; params = [])`
Computes the actual firm payoffs.  Uses parameters computed from one run of the LTE.
"""
function SinglePay(s::simh,
                    mpats::ProjectModule.patientcount,
                    ppats::ProjectModule.patientcount;
                    alf1::Float64 = 829.49,
                    alf2::Float64 = 36166.6,
                    alf3::Float64 = 16309.47,
                    gamma_1_385::Float64 = 20680.0, # ✓
                    gamma_2_385::Float64 = 42692.37, # ✓
                    gamma_3_385::Float64 = 20962.97, # ✓
                    gamma_1_386::Float64 = 81918.29, # X
                    gamma_2_386::Float64 = 74193.4, # X
                    gamma_3_386::Float64 = 99065.79, # X
                    gamma_1_387::Float64 = 30405.32, # X
                    gamma_2_387::Float64 = 49801.84, # X
                    gamma_3_387::Float64 = 22376.8, # X
                    gamma_1_388::Float64 = 10051.55, # ✓
                    gamma_2_388::Float64 = 19019.18, # X
                    gamma_3_388::Float64 = 33963.5, # X
                    gamma_1_389::Float64 = 29122.89, # X
                    gamma_2_389::Float64 = 14279.58, # X
                    gamma_3_389::Float64 = 20708.15, # X
                    gamma_1_390::Float64 = 22830.05, # X
                    gamma_2_390::Float64 = 6754.76, # X
                    gamma_3_390::Float64 = 3667.42, # ✓
                    gamma_1_391::Float64 = 9089.77, # X
                    gamma_2_391::Float64 = 8120.85, # X
                    gamma_3_391::Float64 = 1900.5, # ✓
                    mcaid385::Float64 = 151380.0,
                    mcaid386::Float64 = 48417.0,
                    mcaid387::Float64 = 18845.0,
                    mcaid388::Float64 = 7507.0,
                    mcaid389::Float64 = 9424.0,
                    mcaid390::Float64 = 4623.0,
                    mcaid391::Float64 = 3664.0) # to DRG mean added 3094 - avg reimbursement for DRGs 370-375 under TX Medicaid (2012)
    outp::Float64 = 0.0
    wtp::Float64 = FindWTP(s)
    if s.level == 1
      outp = alf1*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_1_385*(ppats.count385+mpats.count385) - gamma_1_386*(ppats.count386+mpats.count386) - gamma_1_387*(ppats.count387+mpats.count387) - gamma_1_388*(mpats.count388+ppats.count388) - gamma_1_389*(mpats.count389+ppats.count389) - gamma_1_390*(ppats.count390+mpats.count390) - gamma_1_391*(ppats.count391+mpats.count391)
    elseif s.level == 2
      outp = alf2*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_2_385*(ppats.count385+mpats.count385) - gamma_2_386*(ppats.count386+mpats.count386) - gamma_2_387*(ppats.count387+mpats.count387) - gamma_2_388*(mpats.count388+ppats.count388) - gamma_2_389*(mpats.count389+ppats.count389) - gamma_2_390*(ppats.count390+mpats.count390) - gamma_2_391*(ppats.count391+mpats.count391)
    else # level is 3
      outp = alf3*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_3_385*(ppats.count385+mpats.count385) - gamma_3_386*(ppats.count386+mpats.count386) - gamma_3_387*(ppats.count387+mpats.count387) - gamma_3_388*(mpats.count388+ppats.count388) - gamma_3_389*(mpats.count389+ppats.count389) - gamma_3_390*(ppats.count390+mpats.count390) - gamma_3_391*(ppats.count391+mpats.count391)
    end
    return outp
end



"""
`StartingVals(h::simh)`s
This will generate a set of starting values for the hospitals in the market.
Idea - just take the PDV of one period's profit.  But it needs to be done
over the whole set of possible actions.  Right now this assumes that *all* of the
actions have the same value, except for exit.  That can be improved.
"""
function StartingVals(h::simh,
                      ppats::patientcount,
                      mpats::patientcount;
                      disc::Float64 = 0.95)
  return vcat(repmat([max(SinglePay(h, ppats, mpats)/(1-disc), 100000.0)],3), [max(SinglePay(h, ppats, mpats)/((1-disc)*1000), 1000.0)])
end




"""
`ComputeR(hosp::simh, ppats::Dict{Int64, ProjectModule.patientcount}, mpats::Dict{Int64, ProjectModule.patientcount}, action::Int64, iterations::Int64;disc::Float64 = 0.95 )`
Computes the return (current profit + expected continuation) for each hospital in the state.
"""
function ComputeR(hosp::simh,
                  ppats::patientcount,
                  mpats::patientcount,
                  action::Int64,
                  iterations::Int64;
                  disc::Float64 = 0.95)
  try hosp.visited[KeyCreate(hosp.cns, hosp.level)]
    wt::Float64 = 1.0
    if iterations <= 20_000_000
      wt = 1/sqrt(hosp.visited[KeyCreate(hosp.cns, hosp.level)].counter[action])
    else
      wt = 1/hosp.visited[KeyCreate(hosp.cns, hosp.level)].counter[action]
    end
    hosp.visited[KeyCreate(hosp.cns, hosp.level)].aw[action] = (wt)*(SinglePay(hosp, ppats, mpats) + disc*(WProb(hosp.visited[KeyCreate(hosp.cns, hosp.level)]))) + (1-wt)*(hosp.visited[KeyCreate(hosp.cns, hosp.level)].aw[action])
    for el in 1:size(hosp.visited[KeyCreate(hosp.cns, hosp.level)].psi[1,:],1)
      if hosp.visited[KeyCreate(hosp.cns, hosp.level)].psi[1,:] == action
        hosp.visited[KeyCreate(hosp.cns, hosp.level)].psi = hosp.visited[KeyCreate(hosp.cns, hosp.level)].aw[action] # this will change the value below so that the update changes something.
      end
    end
    hosp.visited[KeyCreate(hosp.cns, hosp.level)].psi[2,:] = DA(hosp.visited[hosp.cns].aw)
    hosp.visited[KeyCreate(hosp.cns, hosp.level)].counter[action] += 1
  catch y
    if isa(y, KeyError)
      hosp.visited[KeyCreate(hosp.cns, hosp.level)]=nlrec(  MD(ChoicesAvailable(hosp), StartingVals(hosp, ppats, mpats))  , vcat(ChoicesAvailable(hosp),transpose(PolicyUpdate(StartingVals(hosp, ppats, mpats)))), Dict(k => 1 for k in ChoicesAvailable(hosp)) )
      hosp.visited[KeyCreate(hosp.cns, hosp.level)].counter[action] += 1
    else
      return y
    end
  end
end

"""
`PolicyUpdate(hosp::simh, neww::Array{Float64,1})`
Takes an array of the new W values and maps back to the simh action probabilities.
"""
function PolicyUpdate(neww::Array{Float64,1})
  return exp(neww - maximum(neww))/sum(exp(neww-maximum(neww)))
end

"""
`MD`
makes a dictionary of a certain kind out of array elements.  Maybe this is
in the language or doable with a comprehension.
"""
function MD(a1::Array{Int64,2}, a2::Array{Float64,1})
  d = Dict{Int64,Float64}()
  for el in 1:size(a1,2)
    d[a1[el]] = a2[el]
  end
  return d
end


"""
`DA(d::Dict{Int64,Float64})`
Takes a dictionary and returns PolicyUpdate applied to the elements
"""
function DA(d::Dict{Int64, Float64})
  names::Array{Int64,1}=Array{Int64,1}()
  outp::Array{Float64,1}=Array{Float64,1}()
  for k in keys(d)
    push!(names, k)
    push!(outp, d[k])
  end
  return transpose(hcat(names,  PolicyUpdate(outp)))
end








function CheckConvergence()


end
