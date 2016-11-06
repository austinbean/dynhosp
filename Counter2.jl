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
-  exit::Bool # did it exit?
-  tbu::Bool # Does the record need to be updated ?
-  converged::Bool # has the hospital converged or not?
#NB: TODO - add a level last period value so I can record when the charges for changing need to be assessed.
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
  converged::Bool
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
TexasEq = MakeNew(ProjectModule.fips, ProjectModule.alldists);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists);
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
                     false,
                     false) # added "Converged" Bool.
      for k2 in keys(Tex.mkts)
        for hk2 in keys(Tex.mkts[k2].collection)
          #TODO - note that this might exclude some neighbors from the same county which I want
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
  for el in outp.all
    # Creates an initial value in the "visited" states container.
    el.visited[KeyCreate(el.cns, el.level)] = nlrec(MD(ChoicesAvailable(el), StartingVals(el, ProjectModule.patientcount(5,6,4,13,8,41,248), ProjectModule.patientcount(5,6,4,13,8,41,248))),
                                                    vcat(ChoicesAvailable(el),transpose(PolicyUpdate(StartingVals(el, ProjectModule.patientcount(5,6,4,13,8,41,248), ProjectModule.patientcount(5,6,4,13,8,41,248))))),
                                                    Dict(k => 1 for k in ChoicesAvailable(el)))
    el.visited[KeyCreate(el.cns, el.level)].counter[10] += 1
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
  for el in s.ns
    if el.level != -999
      action = sample(ChoicesAvailable(el), el.chprobs )                            # Take the action
      di = distance(el.lat, el.long, s.lat, s.long)
      if action == 1
        el.level = 3
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs = WeightVec(vec(logitest((0,1), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 2
        el.level = 2
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((1,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 3
        el.level = 2
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((1,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 4
        el.level = 1
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((0,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 5
        el.level = 1
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((0,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 6
        el.level = 3
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =WeightVec(vec(logitest((0,1), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 11
        el.level = -999
        el.tbu = true
      else # action == 10 - "do nothing" (recall that 7,8,9 were entry actions)
        #nothing.
      end
    end
    FixNN(s)   # this will correct the "neighbors of the neighbors"
    FixMainNs(s)
  end
end

"""
`FixMainNs(s::simh)`
Fix the neighbrs of the main facility.
"""
function FixMainNs(s::simh)
  s.cns = neighbors(0,0,0,0,0,0,0,0,0)
  for el in s.ns
    d1::Float64 = distance(s.lat, s.long, el.lat, el.long)<5
    if (d1<5)
      if el.level == 1
        s.cns.level105 += 1
      elseif el.level == 2
        s.cns.level205 += 1
      elseif el.level == 3
         s.cns.level305 += 1
      else # exited
        #do nothing
      end
    elseif (d1>=5)&(d1<15)
      if el.level == 1
        s.cns.level1515 += 1
      elseif el.level == 2
        s.cns.level2515 += 1
      elseif el.level == 3
         s.cns.level3515 += 1
      else # exited
        #do nothing
      end
    elseif (d1>=15)&(d1<25)
      if el.level == 1
        s.cns.level11525 += 1
      elseif el.level == 2
        s.cns.level21525 += 1
      elseif el.level == 3
         s.cns.level31525 += 1
      else # exited
        #do nothing
      end
    end
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
# NB: TODO - Need to subtract the cost of changing here.
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
    #TODO - the level change tracker will let me know when I need to assess the fixed costs.
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
                  disc::Float64 = 0.95,
                  debug::Bool = true)
  k1::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64} = KeyCreate(hosp.cns, hosp.level)
  if haskey(hosp.visited, k1)
    wt::Float64 = 1.0
    if iterations <= 20_000_000
      wt = 1/sqrt(hosp.visited[k1].counter[action])
    else
      wt = 1/hosp.visited[k1].counter[action]
    end
    hosp.visited[k1].aw[action] = (wt)*(SinglePay(hosp, ppats, mpats) + disc*(WProb(hosp.visited[k1]))) + (1-wt)*(hosp.visited[k1].aw[action])
    hosp.visited[k1].psi = ProbUpdate(hosp.visited[k1].aw)
    hosp.visited[k1].counter[action] += 1
    if debug
      PrintVisited(hosp)
    end
  else # Key not there.
    println("hi")
    hosp.visited[k1]=nlrec(MD(ChoicesAvailable(hosp), StartingVals(hosp, ppats, mpats)), vcat(ChoicesAvailable(hosp),transpose(PolicyUpdate(StartingVals(hosp, ppats, mpats)))), Dict(k => 1 for k in ChoicesAvailable(hosp)) )
    hosp.visited[k1].counter[action] += 1
  end
end


"""
ProbUpdate(aw::Dict{Int64,Float64})
This should update the probabilities.
"""
function ProbUpdate(aw::Dict{Int64,Float64})
  tot::Float64 = 0.0
  outp::Array{Float64,1} = Array{Float64,1}()
  labs::Array{Int64,1} = Array{Int64, 1}()
  for el in keys(aw)
    tot += aw[el]
    push!(outp, aw[el])
    push!(labs, el)
  end
  return transpose(hcat(labs,  PolicyUpdate(outp)   ))
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


"""
`WProb(hosp::simh)`
Compute the return R = π + β ∑ Wᵏ(j,xᵏ) Ψᵏ(j, xᵏ+1 ) + β E [ ϵ | xᵏ+1, Ψᵏ], so this
is the function that will compute the second term: β ∑ Wᵏ(j,xᵏ) Ψᵏ(j, xᵏ+1 ).
"""
function WProb(n::nlrec)
  prd::Float64 = 0.0
  for el in keys(n.aw)
    prd += n.aw[el]*n.psi[2,findfirst(n.psi[1,:], el)]
  end
  return prd
end


"""
`ValApprox(D::DynState)`
This computes the dynamic simulation across all of the facilities in all of the markets.
# TODO - Is GetProb updating the ns firms with the state of simh?
# TODO - probs w/in records are not obviously getting updated by ComputeR.
"""
function ValApprox(D::DynState, itlim::Int64; chunk::Array{Int64,1} = collect(1:size(D.all,1)), debug::Bool = true)
  iterations::Int64 = 0
  converged::Bool = false
  a::ProjectModule.patientcount = patientcount(0,0,0,0,0,0,0)
  b::ProjectModule.patientcount = patientcount(0,0,0,0,0,0,0)
  act::Int64 = 0
  while (iterations<itlim)&&(!converged)
    for el in D.all[chunk]
      if !el.converged                                           # only keep simulating with the ones which haven't converged
        act = ChooseAction(el)                                   # Takes an action and returns it.
        a, b = DSim(el.mk, el.fid)                               # Demand as a result of actions.
        GetProb(el)                                              # action choices by other firms
        ComputeR(el, a, b, act, iterations; debug = debug)
        ExCheck(el)
        if debug
          PrintVisited(el)
        end
      end
    end
    iterations += 1
  end
  converged = Halt(D)
end



"""
`Halt(D::DynState)`
For each firm in the simulation, check whether it has converged or not.  Return false
when one firm has been found which hasn't.
"""
function Halt(D::DynState)
  b::Bool = true
  for el in D.all
    if !el.converged
      b = false
      break
    end
  end
  return b
end

"""
`GetChunk(D::DynState, lim::Int64; lower::Bool = false)`
Takes the DynState type and returns a collection of indices of D.all where the markets are a
particular size or range of sizes less than `lim`.  When `lower` is false, takes the `lim` as a
lower bound and then returns everything with a greater market size.
Returns: `Array{Int64, 1}`
"""
function GetChunk(D::DynState, lim::Int64; lower::Bool = false)
  outp::Array{Int64, 1} = Array{Int64,1}()
  for el in 1:size(D.all,1)
    if !lower
      if sum(D.all[el].cns)<=lim
        push!(outp, el)
      end
    else
      if sum(D.all[el].cns)>lim
        push!(outp, el)
      end
    end
  end
  return outp
end



"""
`ChooseAction(h::simh)`
Chooses the action and returns it.
Uses the choice probs implied by the estimated value functions at the state.
"""
function ChooseAction(h::simh)
  act = convert(Int64, sample(h.visited[KeyCreate(h.cns, h.level)].psi[1,:], WeightVec(h.visited[KeyCreate(h.cns, h.level)].psi[2,:])))
  if act == 11
    h.exit = true
  end
  return act
end


"""
`function ExCheck(h::simh)`
Check if the firm exited, then restart.
When the firm exits we also reset all of the neighbors.
"""
function ExCheck(h::simh)
  if h.exit
    h.exit = false
    for n in h.ns
      if n.level == -999
        n.level = n.truelevel
      end
    end
  end
end


"""
`PrintVisited(h::hosp)`
Print the keys of the visited group and the probabilities of the choices, since they
might not be getting updated at this point.
"""
function PrintVisited(h::simh; simple::Bool = true)
  if !simple
    for el in keys(h.visited)
      println(el)
      println(h.visited[el].psi)
      println(h.visited[el].counter)
    end
  else
    println(KeyCreate(h.cns, h.level))
    println(h.visited[KeyCreate(h.cns, h.level)].psi)
  end
end







"""
`CheckConvergence(h::simh; draws::Int64 = 100)`
Check the convergence criterion in Collard-Wexler or Pakes McGuire.
This can be checked every million iterations.  When that happens,
reset the counters for every state in "visited".
"""
function CheckConvergence(h::simh; draws::Int64 = 100, demands::Int64 = 10)
  statecount::Int64 = h.visited.count # how many unique states are there?   NB: this is NOT exactly the right number.  This will be GREATER than those visited in the last million, starting with the second million
  outp::Array{Float64,1} = Array{Float64, 1}() # TODO - should this be a scalar?
  totvisits::Int64 = 0
  demd1::Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}} = DemandGroup(h, demands) # Compute an initial demand set.
  for k in keys(h.visited)
    # TODO - rememember - only draw from that state if visited, which should be when the counter is greater than 2.
    for d = 1:draws
      # Draw the VALUE of the action to get the next state.
      nextact::Int64 = convert(Int64, sample(h.visited[k].psi[1,:], WeightVec(h.visited[k].psi[2,:]))]) # LevelFunction takes Int64 argument in second place.
      #TODO - careful about the levelcheck since the demand is different when we change the level.
      if haskey(h.visited[KeyCreate(h.cns, LevelFunction(h, nextact))]) # check neighbors/level pair
        if LevelFunction(h, nextact) == h.level
            # here can use the demd above
        else
          h.level =

        end
      else # when I haven't been there before, must take the initial value.
        hosp.visited[k1]=nlrec(MD(ChoicesAvailable(hosp), StartingVals(hosp, ppats, mpats)), vcat(ChoicesAvailable(hosp),transpose(PolicyUpdate(StartingVals(hosp, ppats, mpats)))), Dict(k => 1 for k in ChoicesAvailable(hosp)) )
      end
    end
  end
  # DONE - Resets the counter - this will keep track of who has been visited in the last 1_000_000 iterations.
  for k1 in keys(h.visited)
    for k2 in keys(h.visited[k1].counter)
      h.visited[k1].counter[k2] = 1 # Reset all counter keys to 1 to track which states visited in last million iterations.
    end
  end
  return outp
end




"""
`DemandGroup(h::simh, totl::Int64)`
To speed the demand computation, just compute a group of `totl` demands
and sample uniformly from them.   This will need to be done when the state
changes.  Perhaps only the state of the main firm.   `psamps` is the
private patients (first element of return) `msamps` is the Medicaid patients,
the second element of the return type.
"""
function DemandGroup(h::simh, totl::Int64)
  psamps::Array{patientcount,1} = Array{patientcount, 1}()
  msamps::Array{patientcount,1} = Array{patientcount, 1}()
  for i = 1:totl
    a, b = DSim(h.mk, h.fid)
    push!(psamps, a)
    push!(msamps, b)
  end
  return psamps, msamps
end


"""
`EasyDemand(Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}})`
Given an input in the form of a tuple of patientcounts, this chooses one each of the private and
medicaid types and returns them.  The first element is the private patients, second element the
Medicaid patients.  This is thousands of times faster.  Use this method instead.
"""
function EasyDemand(inp::Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}})
  return sample(inp[1]), sample(inp[2])
end


"""
`ThreeDemands(h::simh, totl::Int64)`
Speculative - would it be best to have three possible demand sets with one for each level?
This would:
- Compute demand at some level.
- Update to another.
- Compute at that level.
- Compute at the third level.
- Return a tuple of demands for medicaid and private patients.
"""
function ThreeDemands(h::simh, totl::Int64)
  outp1::Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}} = (Array{patientcount,1}(), Array{patientcount,1}())
  outp2::Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}} = (Array{patientcount,1}(), Array{patientcount,1}())
  outp3::Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}} = (Array{patientcount,1}(), Array{patientcount,1}())
  if h.level == 1
    outp1 = DemandGroup(h, totl)
    h.level = 2
    h.tbu = true # NB: setting this to true means that the call to UpdateDUtil will update the utility value, the reset this
    UpdateDUtil(h)
    outp2 = DemandGroup(h, totl)
    h.level = 3
    h.tbu = true
    UpdateDUtil(h)
    outp3 = DemandGroup(h, totl)
    h.level = 1
  elseif h.level == 2
    outp2 = DemandGroup(h, totl)
    h.level = 1
    h.tbu = true
    UpdateDUtil(h)
    outp1 = DemandGroup(h, totl)
    h.level = 3
    h.tbu = true
    UpdateDUtil(h)
    outp3 = DemandGroup(h, totl)
    h.level = 2 # reset level
  elseif h.level == 3
    outp3 = DemandGroup(h, totl)
    h.level = 1
    h.tbu = true
    UpdateDUtil(h)
    outp1 = DemandGroup(h, totl)
    h.level = 2
    h.tbu = true
    UpdateDUtil(h)
    outp2 = DemandGroup(h, totl)
    h.level = 3
  end
  return outp1, outp2, outp3
end




"""
`AvgAggState(h::simh)`
To make the computation faster, compute the average state in the rest of the market.  Draw actions for
all of the neighbors, compute the aggregate state, then take the mean over such states.
"""
function AvgAggState(h::simh, drws::Int64)
  tem::Array{neighbors, 1} = Array{neighbors, 1}()
  interim::neighbors = neighbors(0,0,0,0,0,0,0,0,0)
  for i = 1:drws
    GetProb(h)
    interim = sum(interim, h.cns)
  end
  return neighbors( round(Int64, interim.level105/drws), round(Int64,interim.level205/drws ), round(Int64, interim.level305/drws), round(Int64, interim.level1515/drws), round(Int64, interim.level2515/drws), round(Int64, interim.level3515/drws), round(Int64,interim.level11525/drws ), round(Int64, interim.level21525/drws), round(Int64, interim.level31525/drws) )
end








#=


=#
