# Now solve the dynamic model itself.


# Can put these into a market and then into a state.

# use the type neighbors !



"""
`KeyCreate(n::neighbors, l::Int64)`
takes neighbors and level and returns a tuple - this is a better way (hopefully) to make the keys for the nlrec.
"""
function KeyCreate(n::neighbors, l::Int64)
  return (n.level105, n.level205, n.level305, n.level1515, n.level2515, n.level3515, n.level11525, n.level21525, n.level31525, l)
end


"""
`CounterObjects()`
Makes the counterfactual objects again.
dyn = CounterObjects(10);

"""
function CounterObjects(T::Int64)
  TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, T);
  Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}())
  CMakeIt(Tex, ProjectModule.fips);
  FillState(Tex, ProjectModule.alldists, T);
  patients = NewPatients(Tex);
  return DynStateCreate(TexasEq, Tex, patients);
end



"""
`DynStateCreate(Tex::EntireState, Tex2::EntireState, p::patientcollection )`
Create the dynamic records from the existing state, don't bother doing it from scratch.
And these don't have to be organized by zip.
Use the EntireState from the equilibrium simulation, not the first counterfactual.
Make using
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);

Note that the function takes TWO EntireState arguments.  This is super dumb, but
only one of them (containing hospital types) has the bed counts.
"""
function DynStateCreate( Tex::EntireState, Tex2::EntireState, p::patientcollection )
  outp = DynState(Array{simh,1}())
  for k1 in keys(Tex.mkts)
    for hk in keys(Tex.mkts[k1].collection)
      newsimh = simh(Tex.mkts[k1].collection[hk].fid,
                     Tex.mkts[k1].collection[hk].lat,
                     Tex.mkts[k1].collection[hk].long,
                     Tex.mkts[k1].collection[hk].level,
                     Tex.mkts[k1].collection[hk].level,
                     Tex.mkts[k1].collection[hk].level,
                     convert(Int64, Tex2.mkts[k1].collection[hk].bedcount),
                     neighbors(0,0,0,0,0,0,0,0,0),
                     Array{Int64,1}(),
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
              # FIXME - there is a problem here in that not all neighbors get pushed to this list.  
              push!(newsimh.nfids, Tex.mkts[k2].collection[hk2].fid) # should add neighboring facs to the nfids list
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
                                                    Dict(k => 0 for k in ChoicesAvailable(el)))
    el.visited[KeyCreate(el.cns, el.level)].counter[10] += 1
  end
  return outp
end


"""
`DynPatients(p::patientcollection, f::fid)`
This will take a collection of patients and create a cmkt, which is a vector of
cpats.  That is, it will take every zip for which `f` is an option, then create
the collection of patients for those zips.  Note that `f` is a FID for a hospital.

Note that the patients collection is created with NewPatients which is defined in DataStructs.
#TODO - here is the problem with patient linking.  This is not outputting anything?
Some of these are coming up empty.  About 473/1938.  Check how many zips there are in the data.
Yes, this is a problem.  There are 1742 unique zips in the inpatient data.  Find out what happened here.
To test:

Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists);
patients = NewPatients(Tex);
DynPatients(patients, 4530190);

"""
function DynPatients(p::patientcollection, f::Int64 )
  outp::cmkt = cmkt(f, Array{cpats,1}())
  zpc = PatientFind(p, f) # finds the zip codes
  for el in zpc
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
function DetUtils(z::zipcode; switch::Bool = false)
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
  for i =1:maximum(size(ar))
    if ar[i]!=0.0
      denom += (ar[i] = exp(ar[i]))
    end
  end
  return log.(map(x->(1/(1-x)), ar./denom))
end



"""
`DSim(c::cmkt, f::Int64)`
This is going to take the collection of patients and the fid and figure out
how many people choose it.  It's ok fast, but not really fast.

#NB: this demand is dumb.  And incorrect.  The correct thing to do is for the
firm to take the expectation.  But we know what that looks like.  I just need to
update the deterministic component, compute the probability, and then multiply by the
number of people in the zip.  This needs fixing, but will be much faster.
"""
function DSim(c::cmkt, f::Int64; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ))
  pcount::patientcount = patientcount(0,0,0,0,0,0,0)
  mcount::patientcount = patientcount(0,0,0,0,0,0,0)
  for el in c.m
    # NB: here is a parallel opportunity, maybe?  Sum across the zip codes across cores?  Or a threading opportunity?
    # Threading possible *per loop* perhaps possible using the Threads.@threads annotation on each one?
    siz1 = size(el.putils[2,:],1) #siz1 and siz2 should always be the same.
    siz2 = size(el.mutils[2,:],1)
    temparr = zeros(siz1)
    for i = 1:el.pcounts.count385
      if el.putils[1,indmax( el.putils[2,:] +rand!(d, temparr) )] == f
        pcount.count385 += 1
      end
    end
    for i = 1:el.pcounts.count386
      if el.putils[1,indmax( el.putils[2,:] +rand!(d, temparr) )] == f
        pcount.count386 += 1
      end
    end
    for i = 1:el.pcounts.count387
      if el.putils[1,indmax( el.putils[2,:] +rand!(d, temparr) )] == f
        pcount.count387 += 1
      end
    end
    for i = 1:el.pcounts.count388
      if el.putils[1,indmax( el.putils[2,:] +rand!(d, temparr) )] == f
        pcount.count388 += 1
      end
    end
    for i = 1:el.pcounts.count389
      if el.putils[1,indmax( el.putils[2,:] +rand!(d, temparr) )] == f
        pcount.count389 += 1
      end
    end
    for i = 1:el.pcounts.count390
      if el.putils[1,indmax( el.putils[2,:] +rand!(d, temparr) )] == f
        pcount.count390 += 1
      end
    end
    for i = 1:el.pcounts.count391
      if el.putils[1,indmax( el.putils[2,:] +rand!(d, temparr) )] == f
        pcount.count391 += 1
      end
    end
    #NB: Medicaid patients here:
    for i = 1:el.mcounts.count385
      if el.mutils[1,indmax( el.mutils[2,:] +rand!(d, temparr) )] == f
        mcount.count385 += 1
      end
    end
    for i = 1:el.mcounts.count386
      if el.mutils[1,indmax( el.mutils[2,:] +rand!(d, temparr) )] == f
        mcount.count386 += 1
      end
    end
    for i = 1:el.mcounts.count387
      if el.mutils[1,indmax( el.mutils[2,:] +rand!(d, temparr) )] == f
        mcount.count387 += 1
      end
    end
    for i = 1:el.mcounts.count388
      if el.mutils[1,indmax( el.mutils[2,:] +rand!(d, temparr) )] == f
        mcount.count388 += 1
      end
    end
    for i = 1:el.mcounts.count389
      if el.mutils[1,indmax( el.mutils[2,:] +rand!(d, temparr) )] == f
        mcount.count389 += 1
      end
    end
    for i = 1:el.mcounts.count390
      if el.mutils[1,indmax( el.mutils[2,:] +rand!(d, temparr) )] == f
        mcount.count390 += 1
      end
    end
    for i = 1:el.mcounts.count391
      if el.mutils[1,indmax( el.mutils[2,:] +rand!(d, temparr) )] == f
        mcount.count391 += 1
      end
    end
  end
  return pcount, mcount
end



"""
`WTPNew(c::Array{Float64,2}, arr::Array{Float64,2})`
Another attempt at computing WTP.
Should this take a fid argument too?
NB: this requires an array argument.  It does *not* allocate.  The array will have to be reset to zero
every call.
This will also ignore elements in arr which don't affect choices, since it loops over size(c.putils) only.
The size of arr is not that important since all unused values are zero AND they are ignored.
"""
function WTPNew(c::Array{Float64,2}, arr::Array{Float64,2})
  int_sum::Float64 = 0.0
  for el in 1:size(c,2)
    if c[1,el]!=0.0
      arr[1,el] = c[1,el]
      int_sum += (arr[2,el] = exp(c[2,el]))
    end
  end
  for i=1:size(c,2)
    arr[2,i]/=int_sum
  end
end


"""
`DemComp(inparr::Array{Float64,2}, fid::Int64, c::cpats, p_or_m::Bool)`
Takes a fid.  Finds the index in inparr.  Takes the share from inparr.
Multiplies all of patient values in cpat by that number.
Returns the value.

- FIXME Note: Can't compute demand when # of patients is zero! 
- This is an issue of the shares.  

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

temparr = zeros(2, 12) # 12 is max # of facilities. 


DemComp(dyn.all[2].mk.m[1].putils, temparr, p1, dyn.all[2].fid, dyn.all[2].mk.m[1], true )
"""
# FIXME - there is a share issue here.  The problem is in inparr.  Maybe something goes wrong with WTPnew?  In temparr?
function DemComp(inparr::Array{Float64,2}, temparr::Array{Float64,2}, pp::patientcount, fid::Int64, c::cpats, p_or_m::Bool)
  # NB: inparr is a sub-field of c.  inparr is either c.putils or c.mutils.
  index::Int64 = 0
  counter::Int64 = 0
  for i = 1:size(inparr,2)
    if inparr[1,i] == fid
      index = i #reassign
    end
  end
  WTPNew(inparr, temparr) # updates temparr
  if index!=0 # don't look for a facility that isn't there.
    if p_or_m # if true then private
      for j in c.pcounts
        if counter == 0 # FIXME - every line here: add if pp.count385 > 0  Next section too.  
          pp.count385 += temparr[2,index]*j
          counter += 1
        elseif counter == 1
          pp.count386 += temparr[2,index]*j
          counter += 1
        elseif counter == 2
          pp.count387 += temparr[2,index]*j
          counter += 1
        elseif counter == 3
          pp.count388 += temparr[2,index]*j
          counter += 1
        elseif counter == 4
          pp.count389 += temparr[2,index]*j
          counter += 1
        elseif counter == 5
          pp.count390 += temparr[2,index]*j
          counter += 1
        elseif counter == 6
          pp.count391 += temparr[2,index]*j
          counter += 1
        else
          println("eee")
        end
      end
    else
      for j in c.mcounts
        if counter == 0
          pp.count385 += temparr[2,index]*j
          counter += 1
        elseif counter == 1
          pp.count386 += temparr[2,index]*j
          counter += 1
        elseif counter == 2
          pp.count387 += temparr[2,index]*j
          counter += 1
        elseif counter == 3
          pp.count388 += temparr[2,index]*j
          counter += 1
        elseif counter == 4
          pp.count389 += temparr[2,index]*j
          counter += 1
        elseif counter == 5
          pp.count390 += temparr[2,index]*j
          counter += 1
        elseif counter == 6
          pp.count391 += temparr[2,index]*j
          counter += 1
        else
          println("eee")
        end
      end
    end
  end
  ArrayZero(temparr)
end

"""
`ArrayZero(arr::Array{Float64,2})`
This quickly sets the array used in WTPNew back to zero.
"""
function ArrayZero(arr::Array{Float64,2})
  dim1::Int64, dim2::Int64 = size(arr)
  for i = 1:dim1
    for j = 1:dim2
      arr[i,j] = 0.0
    end
  end
end



"""
`DSimNew(c::cmkt, f::Int64, pcount::patientcount, mcount::patientcount; maxh::Int64 = 12)`
This is going to compute a market share at the level of a zip or zip-drg.
That formula exists.  The shares will be proportional to the deterministic
utility component.  We just need all of the utilities in the market and number of
each patient type in the same.
e.g., es.all[1].mk is the cmkt.
then es.all[1].mk.m[i] is the cpats.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
DSimNew(dyn.all[2].mk, dyn.all[2].fid, p1, p2)

"""
function DSimNew(c::cmkt, f::Int64, pcount::patientcount, mcount::patientcount; maxh::Int64 = 12)
  temparr::Array{Float64,2} = zeros(2, maxh)
  for el in c.m
    DemComp(el.putils, temparr, pcount, f, el, true)
    DemComp(el.mutils, temparr, mcount, f, el, false)
  end
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
    elseif sr.level == 3
      return pcoeffs.distance*d+pcoeffs.distsq*(d^2)+pcoeffs.distbed*(sr.beds*d/100)+pcoeffs.closest*(0)+pcoeffs.inten
    else # should cover anyone who exited.
      return 0.0
    end
  else
    if sr.level == 1
      return mcoeffs.distance*d+mcoeffs.distsq*(d^2)+mcoeffs.distbed*(sr.beds*d/100)+mcoeffs.closest*(0)
    elseif sr.level == 2
      return mcoeffs.distance*d+mcoeffs.distsq*(d^2)+mcoeffs.distbed*(sr.beds*d/100)+mcoeffs.closest*(0)+mcoeffs.inter
    elseif sr.level == 3
      return mcoeffs.distance*d+mcoeffs.distsq*(d^2)+mcoeffs.distbed*(sr.beds*d/100)+mcoeffs.closest*(0)+mcoeffs.inten
    else # covers exits.
      return 0.0
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
      else # action == 10 - "do nothing" (7,8,9 were entry actions - currently unused)
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
    d1::Float64 = distance(s.lat, s.long, el.lat, el.long)
    if (d1<5)&&(d1>0)
      if el.level == 1
        s.cns.level105 += 1
      elseif el.level == 2
        s.cns.level205 += 1
      elseif el.level == 3
         s.cns.level305 += 1
      else # exited
        #do nothing
      end
    elseif (d1>=5)&&(d1<15)
      if el.level == 1
        s.cns.level1515 += 1
      elseif el.level == 2
        s.cns.level2515 += 1
      elseif el.level == 3
         s.cns.level3515 += 1
      else # exited
        #do nothing
      end
    elseif (d1>=15)&&(d1<25)
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
denoted by f.  f is the hospital FID.
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
  d2 = DynStateCreate(CreateEmpty(fi, da))
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
        #return [0.0,0.0,0.0,0.0]
  return vcat(repmat([min(SinglePay(h, ppats, mpats, 10)/(1-disc), 100.0)],3), [min(SinglePay(h, ppats, mpats, 10)/((1-disc)*1000), 100.0)])
end




"""
ProbUpdate(aw::Dict{Int64,Float64})
This should update the probabilities.
"""
function ProbUpdate(aw::Dict{Int64,Float64})
  outp::Array{Float64,1} = Array{Float64,1}()
  labs::Array{Int64,1} = Array{Int64, 1}()
  for el in keys(aw)
    push!(outp, aw[el])
    push!(labs, el)
  end
  return transpose(hcat(labs,  PolicyUpdate(outp)   ))
end


"""
`PolicyUpdate(neww::Array{Float64,1}; ep::Float64 = 0.0001)`
Takes an array of the new W values and maps back to the simh action probabilities.
NB: When probs are 0, continuation value of error is problematic, since this is
given by the log, so the value is constrained to be this small positive value.
The problem is basically underflow: when returns are really high to staying in and
much smaller to getting out, the estimated prob is zero.


"""
function PolicyUpdate(neww::Array{Float64,1}; ep::Float64 = 0.000000001)
  return max.(exp.(neww-maximum(neww)), ep)/sum(exp.(neww-maximum(neww))) #6devfix
end



"""
`WeightedProbUpdate(aw::Dict{Int64,Float64}, its::Int64)`
There is an underflow problem in the update of probabilities.  Weight them according to the
number of iterations in the update.
#NB: NOT IN USE at the moment.
"""
function WeightedProbUpdate(aw::Dict{Int64, Float64}, ps::Array{Float64,2}, its::Int64)
  next::Array{Float64,2} = ProbUpdate(aw)
  for el1 in 1:maximum(size(ps[1,:]))
    for el2 in 1:maximum(size(next[1,:]))
      if ps[1,el1] == next[1,el2]
        ps[2,el1] = (1/its)*(ps[2,el1]) +(1-(1/its))*(next[2,el1])
      end
    end
  end
  return ps
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
`WProb(n::nlrec)`
Compute the return R = π + β ∑ Wᵏ(j,xᵏ) Ψᵏ(j, xᵏ+1 ) + β E [ ϵ | xᵏ+1, Ψᵏ], so this
is the function that will compute the second term: β ∑ Wᵏ(j,xᵏ) Ψᵏ(j, xᵏ+1 ).
This is the continuation value ignoring the error.
"""
function WProb(n::nlrec)
  prd::Float64 = 0.0
  for el in keys(n.aw)
    prd += n.aw[el]*n.psi[2,findfirst(n.psi[1,:], el)]
  end
  return prd
end



"""
`ContError(n::nlrec)`
Computes the continuation value of the error.
"""
function ContError(n::nlrec)
  return (eulergamma - dot(log.(n.psi[2,:]), n.psi[2,:])) #6devfix
end


"""
`WhyNaN(n::nlrec)`
Getting some NaNs in both probabilities and values in nlrecs.
Try to diagnose that with this.
"""
function WhyNaN(n::nlrec)
  for el in keys(n.aw)
    if isnan(n.aw[el])
      println(n.aw[el], "  ", el)
    end
  end
  for el in n.psi[2,:]
    if isnan(el)
      println(el)
    end
  end
end



"""
`SinglePay(s::simh, mpats::ProjectModule.patientcount, ppats::ProjectModule.patientcount; params = [])`
Computes the actual firm payoffs.  Uses parameters computed from one run of the LTE.
"""
function SinglePay(s::simh,
                    mpats::ProjectModule.patientcount,
                    ppats::ProjectModule.patientcount,
                    action::Int64;
                    scalefact::Float64 = 3.0e9,
                    alf1::Float64 = 8336.17,
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
                    level12::Float64 = 1.64669492e6,
                    level13::Float64 = 5.0165876e6,
                    level21::Float64 = -366430.33,
                    level23::Float64 = 1.83969306e6,
                    level31::Float64 = -90614.32,
                    level32::Float64 = -157206.98,
                    mcaid385::Float64 = 151380.0,
                    mcaid386::Float64 = 48417.0,
                    mcaid387::Float64 = 18845.0,
                    mcaid388::Float64 = 7507.0,
                    mcaid389::Float64 = 9424.0,
                    mcaid390::Float64 = 4623.0,
                    mcaid391::Float64 = 3664.0) # to DRG mean added 3094 - avg reimbursement for DRGs 370-375 under TX Medicaid (2012)
    outp::Float64 = 0.0
    levelc::Float64 = 0.0
    wtp::Float64 = FindWTP(s)
    if s.level == 1&(action!=11)
      outp = alf1*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_1_385*(ppats.count385+mpats.count385) - gamma_1_386*(ppats.count386+mpats.count386) - gamma_1_387*(ppats.count387+mpats.count387) - gamma_1_388*(mpats.count388+ppats.count388) - gamma_1_389*(mpats.count389+ppats.count389) - gamma_1_390*(ppats.count390+mpats.count390) - gamma_1_391*(ppats.count391+mpats.count391)
    elseif s.level == 2&(action!=11)
      outp = alf2*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_2_385*(ppats.count385+mpats.count385) - gamma_2_386*(ppats.count386+mpats.count386) - gamma_2_387*(ppats.count387+mpats.count387) - gamma_2_388*(mpats.count388+ppats.count388) - gamma_2_389*(mpats.count389+ppats.count389) - gamma_2_390*(ppats.count390+mpats.count390) - gamma_2_391*(ppats.count391+mpats.count391)
    elseif s.level == 3&(action!=11)
      outp = alf3*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_3_385*(ppats.count385+mpats.count385) - gamma_3_386*(ppats.count386+mpats.count386) - gamma_3_387*(ppats.count387+mpats.count387) - gamma_3_388*(mpats.count388+ppats.count388) - gamma_3_389*(mpats.count389+ppats.count389) - gamma_3_390*(ppats.count390+mpats.count390) - gamma_3_391*(ppats.count391+mpats.count391)
    else # level = -999 (exit)
      outp = 0.0
    end
    if s.level != s.previous
      if s.level == 1
        if s.previous == 2
          levelc = level21
        elseif s.previous == 3
          levelc = level31
        else
          #not possible
        end
      elseif s.level == 2
        if s.previous == 1
          levelc = level12
        elseif s.previous == 3
          levelc = level32
        else
          # not possible
        end
      elseif s.level == 3
        if s.previous == 2
          levelc = level32
        elseif s.previous == 1
          levelc = level31
        else
          # not possible
        end
      elseif (s.level == -999)||(action == 11) #TODO - here add the scrap values of exiting.
        levelc = 0.0
      end
    end
    return (outp - levelc)/scalefact
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
                  debug::Bool = false)
  k1::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64} = KeyCreate(hosp.cns, hosp.level)
  if haskey(hosp.visited, k1)
    wt::Float64 = 1.0
    if iterations <= 20_000_000
      wt = 1/sqrt(hosp.visited[k1].counter[action]+1)
    else
      wt = 1/(hosp.visited[k1].counter[action]+1)
    end
    hosp.visited[k1].aw[action] = (wt)*(SinglePay(hosp, ppats, mpats, action) + disc*(WProb(hosp.visited[k1])) ) + (1-wt)*(hosp.visited[k1].aw[action])
    if action != 11 # there is no continuation value of the error when the firm exits.
      hosp.visited[k1].aw[action] += disc*ContError(hosp.visited[k1])
    end
    hosp.visited[k1].psi = ProbUpdate(hosp.visited[k1].aw) #WeightedProbUpdate(hosp.visited[k1].aw, hosp.visited[k1].psi, iterations)
    hosp.visited[k1].counter[action] += 1
    hosp.previous = hosp.level # need to record when the level changes.
  else # Key not there.
    if (hosp.level != -999)&&(action != 11) # nothing added for exiters.
      hosp.visited[k1]=nlrec(MD(ChoicesAvailable(hosp), StartingVals(hosp, ppats, mpats)), vcat(ChoicesAvailable(hosp),transpose(PolicyUpdate(StartingVals(hosp, ppats, mpats)))), Dict(k => 0 for k in ChoicesAvailable(hosp)) )
      hosp.visited[k1].counter[action] += 1
    end
  end
end



"""
`RTuple(h::simh, a::Int64)`
Makes a tuple out of a record, level and action.
"""
function RTuple(h::simh, a::Int64)
  return t::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64} = (h.cns.level105, h.cns.level205, h.cns.level305, h.cns.level1515, h.cns.level2515, h.cns.level3515, h.cns.level11525, h.cns.level21525, h.cns.level31525, h.level, a)
end


"""
`PatientZero(mc::patientcount, pc::patientcount)`
Simple function - returns the values of the patientcounts to zeros.
Faster to hard code this for exactly two patientcounts than to call it twice.
"""
function PatientZero(mc::patientcount,  pc::patientcount)
  mc.count385 = 0.0
  mc.count386 = 0.0
  mc.count387 = 0.0
  mc.count388 = 0.0
  mc.count389 = 0.0
  mc.count390 = 0.0
  mc.count391 = 0.0
  pc.count385 = 0.0
  pc.count386 = 0.0
  pc.count387 = 0.0
  pc.count388 = 0.0
  pc.count389 = 0.0
  pc.count390 = 0.0
  pc.count391 = 0.0
end



#### Below this line... Exact Value development ### 


"""
`ExactVal(D::DynState, V::allvisits, itlim::Int64, chunk::Array{Int64,1}; debug::Bool = true)`
Computes the exact solution for smaller markets.  1 - 5 firms at most.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);

entries to consideR: 1-15 are all.
Here number of neighbors and entry (dyn.all[x])
Duopoly:
1
4
5
9

0 6 (monopoly)
0 7 (monopoly)
0 8 (monopoly)

triopoly:
18
20
21
24
25
26

4 13
4 14
4 15

#NB: consider an "itlim" ceiling

# testing: 

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients);
ch = [1] # first element
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
ExactVal(dyn, ch, p1, p2)

PatientZero(p1, p2)

NOTES on current problems:
- 04/25: what I need is something like: full update and partial update.  
- some firms  are getting added... that is, fids are getting added which I don't want added.  Where does that happen?  

"""



##### PLACEHOLDER ######




"""
`ExactChoice(temp::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                     stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, 
                     fid::Int64, 
                     location::Int64,
                     p1::patientcount,
                     p2::patientcount,
                     competitors::Array{Int64,1},
                     D::DynState; )`
What action should the firm choose?
Takes two dictionaries, the DynState, computes the best action, returns the value of the action.
Needs to: 
- compute the demand in expectation at EACH possible level.
- compute the profit at EACH possible level.
- state will be recorded in the dyn record. 
- But the key thing is: return the VALUE of the state.   
NB - level won't change.  I can compute the value of being in all of these states depending on the level.
 Ok - the thing is that this must be done for Both facilities and their neighbors, but the notion of the state 
 for neighbors is different. This is important. 


##### TESTING ######
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients); 

# To Run:

d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
#fid = 3490795;
#location = 1;

d1[dyn.all[1].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[1].fid][StateKey(dyn.all[1], dyn.all[1].level)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 2)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 3)] = 0.0


ExactChoice(d1, d2, dyn.all[1].fid, 1, p1, p2,  dyn)
d1[dyn.all[1].fid]


EXTRA: 

                # FIXME - here is a problem.  Keys are being added to these dicts inconsistently.  Don't add all of them.  Why 
                # are these needed anyway? 
                # I want to *not* add these.  What will that break?   
                # if !haskey(stable, fid) # this should not be necessary when this is debugged.  
                #   stable[fid] = Dict{NTuple{10, Int64},  Float64 }()
                # end 
                # for el in keys(recs) # this adds a record for each of the (state,level) options.  They are put in the stable dict.  
                #   if !haskey(stable, el)
                #     stable[el] = Dict{NTuple{10,Int64}, Float64}()
                #   end
                #   if !haskey( stable[el], TAddLevel(recs[el], 1) )
                #     stable[el][TAddLevel(recs[el], 1)] = 0.5
                #   end 
                #   if !haskey( stable[el], TAddLevel(recs[el], 2) )
                #     stable[el][TAddLevel(recs[el], 2)] = 0.5
                #   end 
                #   if !haskey( stable[el], TAddLevel(recs[el], 3) )
                #     stable[el][TAddLevel(recs[el], 3)] = 0.5
                #   end 
                # end 

    # FIXME - I don't want to FindComps when this is a neighbor.
    # should I permanently take as arguments all of the neighbors I want to do?  
    # think about the 2-3 firm case... especially with non-overlapping sets of neighbors.  
"""

function ExactChoice(temp::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                     stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, 
                     nbs::Dict{Int64, Bool},
                     fid::Int64, 
                     location::Int64,
                     p1::patientcount,
                     p2::patientcount,
                     D::DynState; 
                     messages::Bool = true, 
                     β::Float64 = 0.95,
                     ϕ13::Float64 = 0.0, # substitute values and scale these!
                     ϕ12::Float64 = 0.0,
                     ϕ1EX::Float64 = 0.0,
                     ϕ23::Float64 = 0.0,
                     ϕ21::Float64 = 0.0,
                     ϕ2EX::Float64 = 0.0,
                     ϕ31::Float64 = 0.0,
                     ϕ32::Float64 = 0.0,
                     ϕ3EX::Float64 = 0.0)  

    neighbors::Array{Int64,1} = FindComps(D.all[location], D) # find the competitors.  
    # FIXME 04/25/2017 - the call to ContProbs already fs up because its looking for something that isn't there. 
    # but it isn't looking for it in stable.  it's in nlocs. 
    println("Exact Choice called on: ", fid)
    nloc = Array{Int64,1}() # initialize empty.
    if nbs[D.all[location].fid] # this will be true when this is a "neighbor" only, not a real facility.  
      nloc = FindComps(D.all[location], D) 
    else
      nloc = Array{Int64,1}() # I guess this array will be empty when the firm is the neighbor.  
    end 
    # FIXME - there are two calls to this function.  I don't seem to use "recs"
    # Also this is getting neighbors of the neighbor, which I do not want.  
    recs = StateRecord(neighbors, location, D) # this returns the dict of the state for the firm whose value is being computed from point of view of the main firm.
    rec = StateRecord(D.all[location].nfids, location, D)     # this is computing the state from the point of view of... the main fac. 
    println("From Exact Choice ") 
    println("keys of stable: ", keys(stable)) 
    println(" temp keys: ",keys(temp))
# 04/26/2017 FIXME - the extra fid is here.  This searches in ContProbs USING rec.  
    println("state record ", rec)   
    println("nlocs: ", nloc)
    for el in nloc 
      println("fid is: ", D.all[el].fid)
    end 
    println("the fid: ", fid)
    println("set of neighbors: ", neighbors)
    println("Exact Choice nbs: ", nbs)
    cps::Dict{Int64,Array{Float64,1}} = ContProbs(rec, nloc, stable, D)  
    nstates::Dict{NTuple{9,Int64},Float64} = TotalCombine(D, location, D.all[location].nfids, cps)
    CV1::Float64 = ContVal(nstates, fid, stable ,1)
    CV2::Float64 = ContVal(nstates, fid, stable ,2)
    CV3::Float64 = ContVal(nstates, fid, stable ,3)   
    
    if messages 
      
      println("stable keys before", keys(stable))  
      println("CV's: ", CV1, " ", CV2, " ", CV3) 
      println("the max was ", temp[fid][StateKey(D.all[location],1)]) 
      println("the rev was ", PatientRev(D.all[location],p1,p2,10)) 
    end
  # Update value at Level 1
    D.all[location].level = 1
    UpdateD(D.all[location])                                  # updates the utility for a new level 
    DSimNew( D.all[location].mk, fid, p1, p2)                 # Computes the demand for that level. 
    temp[fid][StateKey(D.all[location],1)] = maximum([ϕ1EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([β*(CV1),-ϕ12+β*(CV2),-ϕ13+β*(CV3)])])
    D.all[location].level = D.all[location].actual            # resets the level 
    UtilDown(D.all[location])                                 # resets the utility
    PatientZero(p1, p2)                                       # overwrites the patientcount with zeros 
  # Update value at Level 2 (repeats steps above!)
    D.all[location].level = 2
    UpdateD(D.all[location]) # Updates deterministic part of utility.  
    DSimNew( D.all[location].mk, fid, p1, p2) # Computes the demand.   
    temp[fid][StateKey(D.all[location],2)] = maximum([ϕ2EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ21+β*(CV1),β*(CV2),-ϕ23+β*(CV3)])])
    D.all[location].level = D.all[location].actual
    UtilDown(D.all[location])
    PatientZero(p1, p2)
  # Update value at Level 3
    D.all[location].level = 3
    UpdateD(D.all[location])
    DSimNew( D.all[location].mk, fid, p1, p2) # Computes the demand.
    temp[fid][StateKey(D.all[location],3)] = maximum([ϕ3EX, PatientRev(D.all[location],p1,p2,10)+β*maximum([-ϕ31+β*(CV1),-ϕ32+β*(CV2),β*(CV3)])])
    D.all[location].level = D.all[location].actual
    UtilDown(D.all[location])
    PatientZero(p1, p2)
end 





"""
`ExactConvergence(current::Dict{ Int64, Dict{NTuple{10, Int64}, Dict{Int64, Float64} } }, stable::Dict{ Int64, Dict{NTuple{10, Int64}, Dict{Int64, Float64} } }; toler::Float64 =0.001, debug::Bool = true  )`
This will check convergence.  Does this by measuring the maximum difference at every state/action pair 
for each firm.  Returns a boolean recording convergence, but also returns a list of fids of unconverged facilities.
Operates on two dictionaries: one the permanent ("stable") and the other the temporary ("current")

Testing: 
dyn = CounterObjects(10);
test1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }();
test2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }();
test1[dyn.all[6].fid] = Dict{NTuple{10, Int64},  Float64 }();
test2[dyn.all[6].fid] = Dict{NTuple{10, Int64},  Float64 }();
StateEnumerate(dyn.all[6].cns, test1[dyn.all[6].fid])
StateEnumerate(dyn.all[6].cns, test2[dyn.all[6].fid])

test1[dyn.all[6].fid][(0,0,0,0,0,0,0,0,0,1)] = 20 #assign a value.
totest = Dict{Int64,Bool}()
totest[dyn.all[6].fid] = false 
ExactConvergence(test1, test2, totest)

ExactConvergence(test1, test2, totest; debug = false)
(false, [4450450]) # this is returning "converged" FALSE and the list of the unconverged facilities (in this case only one.)
"""
function ExactConvergence(current::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                          stable::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } },
                          totest::Dict{Int64,Bool}; 
                          messages::Bool = true,
                          toler::Float64 =0.1, 
                          debug::Bool = true )
  if messages println("From Exact Convergence: ") end
  if messages println("test ", keys(totest)) end
  if messages println("current ", keys(current)) end
  if messages println("stable ", keys(stable)) end
  converge::Bool = false 
  diffs::Dict{Int64,Float64} = Dict{Int64,Float64}()   # check only the guys still being done.
  newchunk::Dict{Int64, Bool} = Dict{Int64,Bool}()
  for fid in keys(current)                             # checks a subset ONLY, given by those in "current" whose locations are in chunk.  
    if !totest[fid]                                    # keys in totest for which false (i.e., not converged)
      maxdiff::Float64 = 0.0 
      for state in keys(current[fid])                  # states available to the firm.
        if haskey(stable[fid], state)
          if abs(current[fid][state] - stable[fid][state]) > maxdiff # we want MAX difference.  
            maxdiff = abs(current[fid][state] - stable[fid][state])
          end 
        else 
          if messages println("a state wasn't found ") end # check if this is messing anything up
          if abs(current[fid][state]) > maxdiff        # we want MAX difference.  
            maxdiff = abs(current[fid][state])
            stable[fid][state] = 0.5                   # add a new value at the state if it isn't in the dict.
          end
        end 
      end
      diffs[fid] = maxdiff                             # keep track of the max diff.  
    end 
  end 
  if debug 
    println("current differences ")
    println(diffs)
  end 
  # Check for convergence - look at the maximum difference across states for each firm. 
  if messages println("diffs keys are: ", keys(diffs)) end 
  for k1 in keys(diffs)
    converge = converge&(diffs[k1]<toler)
    if messages println("diffs key is: ", k1) end
    if messages println(keys(newchunk)) end
    if diffs[k1] > toler # not converged yet 
      newchunk[k1] = false 
    else 
      totest[k1] = true 
    end 
  end 
  return converge, newchunk
end 




"""
`ContVal( futures::Dict{NTuple{9,Int64},Float64}, fid::Int64, stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, lev::Int64)`


Takes the states of other firms in the Dict futures, computes continuation vals using the Dict stable, 
then returns a float of the CV.  
futures should be the dict returned by TotalCombine - this is a set of probabilities of various outcomes.
stable is the dict of values - i.e., the continuation values of the firms
level is the level for the firm for which we are computing CV 
fid is the dyn.all[xxx].fid - the firm for which we are computing the CV.



d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

d1[dyn.all[18].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[18].fid][StateKey(dyn.all[18], 1)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 2)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 3)] = 0.0
location2 = FindComps(dyn.all[18], dyn) # locations are 19 and 152
recs2 = StateRecord(dyn.all[18].nfids, 18, dyn)

d1[dyn.all[19].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 1)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 2)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 3)] = 0.0

d1[dyn.all[152].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 1)] = 0.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 2)] = 1.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 3)] = 2.0

cp2 = ContProbs(recs2, location2, d1, dyn)
dyn.all[19].nfids = [672285, 373510] # this correction should not be necessary.
compprobs = TotalCombine(dyn, 18, dyn.all[18].nfids, cp2)

ContVal(compprobs, dyn.all[18].fid, d1, 1)

# FIXME - competitors are being added here who should not be.  

"""

function ContVal(futures::Dict{NTuple{9,Int64},Float64}, 
                 fid::Int64, 
                 stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, 
                 lev::Int64)
  outp::Float64 = 0.0
  if lev == 1 
    for k1 in keys(futures)
      if haskey(stable[fid],TAddLevel(k1,lev) )
        outp += futures[k1]*stable[fid][TAddLevel(k1,lev)]
          # if isnan(futures[k1])
          #  println("futures is nan ", futures[k1])
          # end
          # if isnan(stable[fid][TAddLevel(k1,lev)]) 
          #   println("stable is nan ", stable[fid][TAddLevel(k1,lev)])
          # end 
      else 
        # FIXME - here is the problem.  This is added, but should not be.  See three more times below.  
        stable[fid][TAddLevel(k1,lev)] = 0.5
      end 
    end 
  elseif lev == 2
    for k1 in keys(futures)
      if haskey(stable[fid],TAddLevel(k1,lev) )
        outp += futures[k1]*stable[fid][TAddLevel(k1,lev)]
      else 
        stable[fid][TAddLevel(k1,lev)] = 0.5
      end     end 
  elseif lev == 3
    for k1 in keys(futures)
      if haskey(stable[fid],TAddLevel(k1,lev) )
        outp += futures[k1]*stable[fid][TAddLevel(k1,lev)]
      else 
        stable[fid][TAddLevel(k1,lev)] = 0.5
      end 
    end
  else 
    #do nothing.  
  end 
  return outp 
end 







"""
`ContProbs(state_recs::Dict{Int64,NTuple{9,Int64}},
                  nlocs::Array{Int64,1}, # locations of neighbors.
                  stable_vals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64} },
                  D::DynState)`

Picks out the probabilities of actions from opponents values.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients); 

# To Run:

d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
#fid = 3490795;
#location = 1;

d1[dyn.all[1].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[1].fid][StateKey(dyn.all[1], dyn.all[1].level)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 2)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 3)] = 0.0
location = FindComps(dyn.all[1], dyn)
recs = StateRecord(dyn.all[1].nfids, 1, dyn)

d1[1391330] = Dict{NTuple{10, Int64}, Float64}()
d1[1391330][TAddLevel(recs[1391330], 1)] = 0.0
d1[1391330][TAddLevel(recs[1391330], 2)] = 0.0
d1[1391330][TAddLevel(recs[1391330], 3)] = 0.0

cp = ContProbs(recs, location, d1, dyn)
# should return: Dict{Int64,Array{Float64,1}} with 1 entry: 1391330 => [0.333333, 0.333333, 0.333333]

# need a test for two firms as well.  


d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

d1[dyn.all[18].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[18].fid][StateKey(dyn.all[18], 1)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 2)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 3)] = 0.0
location2 = FindComps(dyn.all[18], dyn) # locations are 19 and 152
recs2 = StateRecord(dyn.all[18].nfids, 18, dyn)

d1[dyn.all[19].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 1)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 2)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 3)] = 0.0

d1[dyn.all[152].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 1)] = 0.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 2)] = 0.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 3)] = 0.0

cp = ContProbs(recs2, location2, d1, dyn)

# should return: Dict{Int64,Array{Float64,1}} with 2 entries:  672285 => [0.333333, 0.333333, 0.333333] 373510 => [0.333333, 0.333333, 0.333333]
# FIXME - this is a problem too.  It's iterating over all of these keys.  Why?  OR - what keys and why do I want them?

"""
function ContProbs(state_recs::Dict{Int64,NTuple{9,Int64}},
                   nlocs::Array{Int64,1}, # locations of neighbors.
                   stable_vals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64} },
                   D::DynState)
  outp::Dict{Int64, Array{Float64,1}} = Dict{Int64, Array{Float64,1}}()
  for el in nlocs # these index the locations of neighbors in the array.
    outp[D.all[el].fid] = exp.([stable_vals[D.all[el].fid][TAddLevel(state_recs[D.all[el].fid],1)], stable_vals[D.all[el].fid][TAddLevel(state_recs[D.all[el].fid],2)], stable_vals[D.all[el].fid][TAddLevel(state_recs[D.all[el].fid],3)]])
    outp[D.all[el].fid]./=(sum(outp[D.all[el].fid]))
    if sum(outp[D.all[el].fid]) == 0
      println("ContProbs problem.")
    end   
  end 
  return outp
end


"""
`PMatch(k1::NTuple{9, Int64}, k2::NTuple{3, Int64}, ind::Int64)`
facilitates partial key match.

k1 = (1,2,3,4,5,6,7,8,9)
PMatch(k1, (1,2,3), 1) # true 
PMatch(k1, (1,2,3), 2) # false 
PMatch(k1, (4,5,6), 3) # false 
"""
function PMatch(k1::NTuple{9, Int64}, k2::NTuple{3, Int64}, ind::Int64)
  if ind == 1
    return (k1[1]==k2[1])&(k1[2]==k2[2])&(k1[3]==k2[3])
  elseif ind == 2
    return (k1[4]==k2[1])&(k1[5]==k2[2])&(k1[6]==k2[3])
  elseif ind == 3
    return (k1[7]==k2[1])&(k1[8]==k2[2])&(k1[9]==k2[3])
  else 
    return false 
  end 
end 



"""
`NewCombine(D::DynState, nlocs::Array{Int64,1},location::Int64, contprobs::Dict{Int64,Array{Float64,1}};actions::Array{Int64,1} = [1, 2, 3])`

Ok. Should take the output of contprobs and return the right output.

d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

d1[dyn.all[18].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[18].fid][StateKey(dyn.all[18], 1)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 2)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 3)] = 0.0
location2 = FindComps(dyn.all[18], dyn) # locations are 19 and 152
recs2 = StateRecord(dyn.all[18].nfids, 18, dyn)

d1[dyn.all[19].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 1)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 2)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 3)] = 0.0

d1[dyn.all[152].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 1)] = 0.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 2)] = 1.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 3)] = 2.0

cp = ContProbs(recs2, location2, d1, dyn)
dyn.all[19].nfids = [672285, 373510] # this correction should not be necessary.  
PartialCombine([672285, 373510], cp )
"""
function PartialCombine(fids::Array{Int64,1},
                        contprobs::Dict{Int64,Array{Float64,1}})
  outp::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  push!(outp, ([0, 0, 0], 1)) #initialize the array.  
  for k1 in fids # locations of competitors
    loc1 = [0.0 0.0 0.0]; loc1[1] = contprobs[k1][1];
    loc2 = [0.0 0.0 0.0]; loc2[2] = contprobs[k1][2];
    loc3 = [0.0 0.0 0.0]; loc3[3] = contprobs[k1][3];
    outp = GenStates(outp, loc1, loc2, loc3)
  end 
  return outp 
end 


"""
`TotalCombine`
This will do what I want.  
Return all of the states, with the associated probabilities.  

d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
#fid = 3490795;
#location = 1;

d1[dyn.all[1].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[1].fid][StateKey(dyn.all[1], dyn.all[1].level)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 2)] = 0.0
d1[dyn.all[1].fid][StateKey(dyn.all[1], 3)] = 0.0
location = FindComps(dyn.all[1], dyn)
recs = StateRecord(dyn.all[1].nfids, 1, dyn)

d1[1391330] = Dict{NTuple{10, Int64}, Float64}()
d1[1391330][TAddLevel(recs[1391330], 1)] = 0.0
d1[1391330][TAddLevel(recs[1391330], 2)] = 0.0
d1[1391330][TAddLevel(recs[1391330], 3)] = 0.0
cp = ContProbs(recs, location, d1, dyn)

TotalCombine(dyn, 1, dyn.all[1].nfids, cp)


# try with two firms...


d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

d1[dyn.all[18].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[18].fid][StateKey(dyn.all[18], 1)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 2)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 3)] = 0.0
location2 = FindComps(dyn.all[18], dyn) # locations are 19 and 152
recs2 = StateRecord(dyn.all[18].nfids, 18, dyn)

d1[dyn.all[19].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 1)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 2)] = 0.0
d1[dyn.all[19].fid][TAddLevel(recs2[dyn.all[19].fid], 3)] = 0.0

d1[dyn.all[152].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 1)] = 0.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 2)] = 1.0
d1[dyn.all[152].fid][TAddLevel(recs2[dyn.all[152].fid], 3)] = 2.0

cp2 = ContProbs(recs2, location2, d1, dyn)
dyn.all[19].nfids = [672285, 373510] # this correction should not be necessary.
compprobs = TotalCombine(dyn, 18, dyn.all[18].nfids, cp2)

"""

function TotalCombine(D::DynState,
                      location::Int64,
                      nfids::Array{Int64,1},
                      contprobs::Dict{Int64,Array{Float64,1}})
  # Create the outputs
  out05::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  push!(out05, ([0, 0, 0], 1.0))
  out515::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  push!(out515, ([0, 0, 0], 1.0))
  out1525::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  push!(out1525, ([0, 0, 0], 1.0))
  lastout::Dict{NTuple{9, Int64}, Float64} = Dict{NTuple{9, Int64}, Float64}()
  # Competitors' locations: 
  comps::Array{Int64,1} = FindComps(D.all[location], D)
  # Measure the distances, apply GenStates to relevant segment.
  loc1 = [0.0 0.0 0.0]; 
  loc2 = [0.0 0.0 0.0];
  loc3 = [0.0 0.0 0.0];   
  for k1 in comps 
    d1::Float64 = distance(D.all[location].lat, D.all[location].long, D.all[k1].lat, D.all[k1].long) # how far from the main fac?
    loc1[1] = contprobs[D.all[k1].fid][1];
    loc2[2] = contprobs[D.all[k1].fid][2];
    loc3[3] = contprobs[D.all[k1].fid][3];
    if (d1>0)&(d1<5)
      out05 = GenStates(out05, loc1, loc2, loc3)
    elseif (d1>=5)&(d1<15)
      out515 = GenStates(out515, loc1, loc2, loc3)
    elseif (d1>=15)&(d1<25)
      out1525 = GenStates(out1525, loc1, loc2, loc3)
    else 
      # do nothing 
    end 
    loc1[1] = 0.0; loc2[2] = 0.0; loc3[3] = 0.0; # reassign to 0.
  end 
  for el1 in out05
    for el2 in out515
      for el3 in out1525
        lastout[Tuple((el1[1]..., el2[1]..., el3[1]...))] = el1[2]*el2[2]*el3[2] 
      end 
    end 
  end 
  return lastout 
end 








"""
`CombineV(args...; len = 3)`
Takes a set of vectors of length 3.
Returns a count of the max element across vectors.
Multiplies that max element and returns as prob.

# Testing:
CombineV( [1/2 0 0], [1/2 0 0], [0 0 2])
# ((2, 0, 1), 0.5)

"""
function CombineV(args...; len = 3)
  outp::Array{Int64,1} = zeros(Int64, 3)
  prob::Float64 = 1.0
  for (i, arg) in enumerate(args)
    if length(arg)==len 
      val, indx = findmax(arg)
      outp[indx] += 1
      prob *= val 
    end 
  end 
  return Tuple(outp), prob 
end 


"""
`CombineVInput(inpt::Array{Int64,1}, pr::Float64, args...)` 
Similar to CombineV but takes a vector and float input 

# testing 
CombineVInput([1, 0, 0], 0.5, [0, 0.3, 0]) == ([1, 1, 0], 0.15)
# 
CombineVInput([1, 0, 0], 0.5, [0, 0.3, 0], [0, 0.3, 0]) == ([1, 2, 0], 0.045)
This modifies the vector inpt. 

#NOTE - this copies the input and returns a new vector.  That may be the best way but does allocate more.   
"""
function CombineVInput(inpt::Array{Int64,1}, pr::Float64, args...)
  outp::Array{Int64,1} = zeros(Int64, size(inpt))
  fl::Float64 = pr
  for i = 1:length(inpt)
    outp[i] = inpt[i]
  end   
  for (i, arg) in enumerate(args)
    if length(arg)==length(inpt)
      val, indx = findmax(arg)
      outp[indx] += 1
      fl *= val 
    end 
  end 
  return outp, fl  
end 


"""
`GenStates(inp::Array{Tuple{Array{Int64,1}, Float64}, 1}, args... )`
Takes a collection of Tuples{ Array{Int64,1}, Float64} and a list of args and 
applies CombineVInput to each element.  Returns a new collection of the same.
The number of elements should be growing with EACH call, according to something like 
 size(inp)*number of "args"


#Testing 

in1 = [([1,0,0], 1.0)]
GenStates(in1, [0,0.5,0], [0,0,0.5], [0.3, 0, 0])

# would like: ([1 1 0], 0.5), ([1 0 1], 0.5), ([2 0 0], 0.3)
# the issue might be that inpt above is getting modified!
# I need to keep that element constant and evaluate these elements on it.  
"""
function GenStates(inp::Array{Tuple{Array{Int64,1}, Float64}, 1}, args... )
  outp::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  for (i, arg) in enumerate(args)
    for el in inp # this is a Tuple, Float 
      push!(outp, CombineVInput(el[1], el[2], arg))
    end 
  end 
  return outp 
end






"""
`StateRecord(neighbors::Array{Int64,1}, location::Int64, D::DynState)`

Ugh.  What this function does is compute the state of neighboring facilities RELATIVE to the organization of the fac for 
which we are computing the eq.  Imagine Hosp A has one neighbor B at 25 miles, but B has more neighbors at 10 further miles 
distant.  I want the state for B to reflect ONLY the existence of A.  I think this is the right thing to do.  
"""
function StateRecord(neighbors::Array{Int64,1}, # an array of fids 
                     location::Int64, # the location of the main hospital in dyn.all[]
                     D::DynState)
  outp::Dict{Int64, NTuple{9,Int64}} = Dict{Int64, NTuple{9,Int64}}() # records the neighbors of competing firms, relative to the firm we are computing EQ for, 
  locs = FindComps(D.all[location], D) # where are the neighbors 
  for el1 in locs # vector of locations of neighbors
      intermed::Array{Int64,1} = zeros(Int64,9)
      d1::Float64 = distance(D.all[location].lat, D.all[location].long, D.all[el1].lat, D.all[el1].long)
      if (d1>0)&(d1<=5)
        if D.all[location].level == 1
          intermed[1]+=1
        elseif D.all[location].level == 2
          intermed[2]+=1
        elseif D.all[location].level == 3
          intermed[3]+=1
        else
          #do nothing 
        end 
      elseif (d1>5)&(d1<=15)
        if D.all[location].level == 1
          intermed[4]+=1
        elseif D.all[location].level == 2
          intermed[5]+=1
        elseif D.all[location].level == 3
          intermed[6]+=1
        else
          #do nothing 
        end 
      elseif (d1>15)&(d1<25)
        if D.all[location].level == 1
          intermed[7]+=1
        elseif D.all[location].level == 2
          intermed[8]+=1
        elseif D.all[location].level == 3
          intermed[9]+=1
        else
          #do nothing 
        end 
      else 
        # do nothing 
      end 
      for el2 in locs 
        if el1!=el2
          d2::Float64 = distance(D.all[el1].lat, D.all[el1].long, D.all[el2].lat, D.all[el2].long)
          if (d2>0)&(d2<=5)
            if D.all[el2].level == 1
              intermed[1]+=1
            elseif D.all[el2].level == 2
              intermed[2]+=1
            elseif D.all[el2].level == 3
              intermed[3]+=1
            else
              #do nothing 
            end 
          elseif (d2>5)&(d2<=15)
            if D.all[el2].level == 1
              intermed[4]+=1
            elseif D.all[el2].level == 2
              intermed[5]+=1
            elseif D.all[el2].level == 3
              intermed[6]+=1
            else
              #do nothing 
            end 
          elseif (d2>15)&(d2<25)
            if D.all[el2].level == 1
              intermed[7]+=1
            elseif D.all[el2].level == 2
              intermed[8]+=1
            elseif D.all[el2].level == 3
              intermed[9]+=1
            else
              #do nothing 
            end 
          else 
            # do nothing 
          end 
        end 
      end
      outp[D.all[el1].fid] = Tuple(intermed)
  end 
  return outp 
end 


"""
`TAddLevel(t1::NTuple{9,Int64}, l::Int64)`
Just adds another element l to the end of the tuple t1.
"""
function TAddLevel(t1::NTuple{9,Int64}, l::Int64)
  return (t1..., l)
end 


"""
`TupletoCNS(NTuple{9, Int64})`
Takes a tuple of Ints and returns the same as type neighbors 
This is for use with the output of StateRecord.
"""
function TupletoCNS(n::NTuple{9,Int64})
  return ProjectModule.neighbors(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9])
end 




"""
`PatientRev(s::simh, mpats::ProjectModule.patientcount, ppats::ProjectModule.patientcount; params = [])`
Computes the actual firm payoffs.  Uses parameters computed from one run of the LTE.
This one does *not* include fixed costs of changing level.  Otherwise it is identical to `SinglePay`

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);


dyn = DynStateCreate(TexasEq, Tex, patients); 

p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

DSimNew(dyn.all[1].mk, dyn.all[1].fid, p1, p2)

PatientRev(dyn.all[1], p1, p2, 10) # recall that the scale factor is large so this should be small.
"""
function PatientRev(s::simh,
                    mpats::ProjectModule.patientcount,
                    ppats::ProjectModule.patientcount,
                    action::Int64;
                    scalefact::Float64 = 3.0e9,
                    alf1::Float64 = 8336.17,
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
    elseif s.level == 3
      outp = alf3*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_3_385*(ppats.count385+mpats.count385) - gamma_3_386*(ppats.count386+mpats.count386) - gamma_3_387*(ppats.count387+mpats.count387) - gamma_3_388*(mpats.count388+ppats.count388) - gamma_3_389*(mpats.count389+ppats.count389) - gamma_3_390*(ppats.count390+mpats.count390) - gamma_3_391*(ppats.count391+mpats.count391)
    else # level = -999 (exit)
      outp = pi*scalefact
    end
    return outp/scalefact
end




"""
`StateKey(h::simh)`
Returns a tuple of the neighbors and level to use in `ExactVal()`
"""
function StateKey(h::simh, a::Int64)
  return (h.cns.level105, h.cns.level205, h.cns.level305, h.cns.level1515, h.cns.level2515, h.cns.level3515, h.cns.level11525, h.cns.level21525, h.cns.level31525, a)
end 


"""
`FindComps(h::simh, D::DynState)`
Finds the locations of the neighbors in the DynState.  
"""
function FindComps(h::simh, D::DynState)
    # takes a state and finds me the neighbors.
    outp::Array{Int64,1} = Array{Int64,1}()
    for el in h.nfids
      for loc in 1:size(D.all,1) 
        if D.all[loc].fid == el 
          push!(outp, loc)
        end  
      end 
    end 
    return outp 
end 


"""
`UpdateD(h::simh)`
Updates the utility component.  The deterministic part.
"""
function UpdateD(h::simh)
    if h.level == 1
      if h.actual == 3
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 3, 1)
        end 
      elseif h.actual == 2
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 2,1) 
        end 
      else # h.actual == 1
        # do nothing.
      end 
    elseif h.level == 2
      if h.actual == 1
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 1, 2)
        end 
      elseif h.actual == 3
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 3, 2)
        end       
      else # h.actual == 2
        # do nothing. 
      end
    else #h.level == 3
      if h.actual == 1
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 1, 3)
        end 
      elseif h.actual == 2
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 2,3)
        end 
      else # h.actual == 3
        # do nothing.
      end
    end 
end 



"""
`UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64;
                 inteninter_med::Float64 = -0.572397,
                 interinten_med::Float64 = 0.572397,
                 inten_med::Float64 = 1.34994,
                 inter_med::Float64 = 0.777542,
                 inteninter_p::Float64 = 0.3197218,
                 interinten_p::Float64 = -0.3197218,
                 inten_p::Float64 = 1.18599,
                 inter_p::Float64 = 0.866268) `
Updates the deterministic component of the utility quickly.  
"""
function UtilUp( c::cpats, 
                 fid::Int64, 
                 actual::Int64, 
                 current::Int64; 
                 inteninter_med::Float64 = -0.572397,
                 interinten_med::Float64 = 0.572397,
                 inten_med::Float64 = 1.34994,
                 inter_med::Float64 = 0.777542,
                 inteninter_p::Float64 = 0.3197218,
                 interinten_p::Float64 = -0.3197218,
                 inten_p::Float64 = 1.18599,
                 inter_p::Float64 = 0.866268)
  indx_m::Int64 = findfirst(c.mutils[1,:], fid)
  indx_p::Int64 = findfirst(c.putils[1,:], fid)
  if (actual == 1)&(current == 1)
    # do nothing.
  elseif (actual == 1)&(current == 2)
    c.putils[2,indx_p] += inter_p
    c.mutils[2,indx_m] += inter_med
  elseif (actual == 1)&(current == 3)
    c.putils[2,indx_p] += inten_p
    c.mutils[2,indx_m] += inten_med 
  elseif (actual == 2)&(current == 1)
    c.putils[2,indx_p] -= inter_p
    c.mutils[2,indx_m] -= inter_med
  elseif (actual == 2)&(current == 2)
    # do nothing.
  elseif (actual == 2)&(current == 3)
    c.putils[2,indx_p] += inteninter_p 
    c.mutils[2,indx_m] += inteninter_med
  elseif (actual == 3)&(current == 1)
    c.putils[2,indx_p] -= inten_p 
    c.mutils[2,indx_m] -= inten_med 
  elseif (actual == 3)&(current == 2)
    c.putils[2,indx_p] -= interinten_p
    c.mutils[2,indx_m] -= interinten_med 
  elseif (actual == 3)&(current == 3)
    # do nothing.
  end 
end 


"""
`UtilDown(h::simh)`
Adjusts the deterministic utility back down.  
"""
function UtilDown(h::simh)
    if h.level == 1
      if h.actual == 3
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 1, 3)
        end 
      elseif h.actual == 2
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 1,2) 
        end 
      else # h.actual == 1
        # do nothing.
      end 
    elseif h.level == 2
      if h.actual == 1
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 2, 1)
        end 
      elseif h.actual == 3
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 2, 3)
        end       
      else # h.actual == 2
        # do nothing. 
      end
    else #h.level == 3
      if h.actual == 1
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 3, 1)
        end 
      elseif h.actual == 2
        for el in 1:size(h.mk.m,1) # this is an array of cpats
          UtilUp(h.mk.m[el], h.fid, 3,2)
        end 
      else # h.actual == 3
        # do nothing.
      end
    end 
end 






"""
`DictClean`
In ExactVal there is a dict which stores the values of current continuations to return 
and another in a temporary.  This cleans the temporary.  

Testing: 
dyn = CounterObjects(10);
test1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }();
test1[dyn.all[6].fid] = Dict{NTuple{10, Int64},  Float64 }();
StateEnumerate(dyn.all[6].cns, test1[dyn.all[6].fid])
test1[dyn.all[6].fid][(0,0,0,0,0,0,0,0,0,1)] = 20 #assign a value.
DictClean(test1)

# should return: 
test1[4450450][(0,0,0,0,0,0,0,0,0,1)] == 0.0

"""
function DictClean(d::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } })
  for k1 in keys(d) # these are fids 
    for k2 in keys(d[k1]) # these are neighbor state/level keys at the hospital level.  
      d[k1][k2] = 0.0 # this should return all values to zero.  
    end 
  end 
end 


"""
`DictCopy`
Copy results from the temporary to the permanent.
Let d1 be the permanent and d2 be the temporary.  

Testing: 
dyn = CounterObjects(10);
test1 = Dict{ Int64, Dict{NTuple{10, Int64},  Float64}  }();
test2 = Dict{ Int64, Dict{NTuple{10, Int64},  Float64}  }();
test1[dyn.all[6].fid] = Dict{NTuple{10, Int64}, Float64 }();
test2[dyn.all[6].fid] = Dict{NTuple{10, Int64},  Float64 }();
StateEnumerate(dyn.all[6].cns, test1[dyn.all[6].fid])
StateEnumerate(dyn.all[6].cns, test2[dyn.all[6].fid])

test1[dyn.all[6].fid][(0,0,0,0,0,0,0,0,0,1)] = 20 #assign a value.

DictCopy(test2, test1)

test2[dyn.all[6].fid][(0,0,0,0,0,0,0,0,0,1)] == 20.0 # should return true

"""
function DictCopy(d1::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, 
                  d2::Dict{ Int64, Dict{NTuple{10, Int64},  Float64}  })
  for k1 in keys(d1) # these are fids 
    for k2 in keys(d1[k1]) # these are neighbor state/level keys at the hospital level.  
      if !haskey(d2[k1],k2 ) # if the temporary doesn't have the key 
        # do nothing???  
      else 
        d1[k1][k2] = d2[k1][k2] # this should copy from the temp to the permanent.   
      end  
    end 
  end 
end 



"""
`StateEnumerate(c::ProjectModule.neighbors, inp::Dict{NTuple{10, Int64},  Float64 }; fxd::Bool = false)`
Enumerates all possible neighbor states.
Note that the number of combinations at each x-x+n distance is given by
the combinatoric ( hospitals + levels - 1, hospitals).  (Google "Stars and Bars")
The total number of state elements is given below by EnumUp(n05) * EnumUp(n515) * EnumUp(n1525) * 3
If fxd == true, then EnumUp( ; fixed = true)

dyn = CounterObjects(10);
outp = Dict{NTuple{10,Int64}, Dict{Int64,Float64}}();
outp[(0,0,0,0,0,0,0,0,0,27)] = Dict(1=>0.0, 2=>0.0); # tests that current entries are not overwritten.
StateEnumerate(dyn.all[1].cns, outp);

# second test:
outp = Dict{NTuple{10,Int64}, Dict{Int64,Float64}}();
outp[(0,0,0,0,0,0,0,0,0,27)] = Dict(1=>0.0, 2=>0.0); # tests that current entries are not overwritten.
StateEnumerate(dyn.all[1].cns, outp; fxd = true);


Will not do exact computations for very large markets, but...
For a big market, see dyn.all[50].cns, which generates around 2 million states.
#TODO - for the future: better guesses here about the initial values, not zero.
"""
function StateEnumerate(c::ProjectModule.neighbors,  inp::Dict{NTuple{10, Int64},  Float64 }; fxd::Bool = false)
  n05::Int64 = c.level105 + c.level205 + c.level305
  n515::Int64 = c.level1515 + c.level2515 + c.level3515
  n1525::Int64 = c.level11525 + c.level21525 + c.level31525
  outp::Dict{NTuple{10, Int64}, Float64 } = Dict{NTuple{10, Int64}, Float64}()
  for i in EnumUp(n05; fixed = fxd)
    for j in EnumUp(n515; fixed = fxd)
      for k in EnumUp(n1525; fixed = fxd)
        if !haskey(outp, TupleSmash(i,j,k,1))
          outp[TupleSmash(i,j,k,1)] = 0.5 #Dict{Int64, Float64}(10=>0.0, 1=>0.0, 2=>0.0, 11=>0.0)
        end
        if !haskey(outp, TupleSmash(i,j,k,2))
          outp[TupleSmash(i,j,k,2)] = 0.5#Dict{Int64, Float64}(5=>0.0, 10=>0.0, 6=>0.0, 11=>0.0)
        end
        if !haskey(outp, TupleSmash(i,j,k,3))
          outp[TupleSmash(i,j,k,3)] = 0.5 #Dict{Int64, Float64}(4=>0.0, 3=>0.0, 10=>0.0, 11=>0.0)
        end
      end
    end
  end
  merge!(inp, outp) # merge the dicts.  Does not return anything explicitly.
end

"""
`EnumerLevel(n::Int64)`
Make the relevant combinations.
Why shouldn't this directly take some tuple arguments?  I know that's what I want
in the end.
"""
function EnumerLevel(n::NTuple{3,Int64})
  firstnonzero::Int64 = findfirst(n)
  t1::NTuple{3,Int64}=(n[1], n[2], n[3]+1)
  t2::NTuple{3,Int64}=(n[1], n[2]+1, n[3])
  t3::NTuple{3,Int64}=(n[1]+1, n[2], n[3])
  if firstnonzero == 0 || firstnonzero==3
    return t1, t2, t3
  elseif firstnonzero == 1
    return (t3,) # NB: must return a tuple else the iteration in StateEnumerate won't work.
  elseif firstnonzero == 2
    return t3,t2
  end
end


"""
`TupleSmash(n1::Tuple{3,Int64}, n2::NTuple{3,Int64}, n3::NTuple{3,Int64}, level::Int64)`
Takes three tuples of ints and combines them into one giant tuple, adds a level at the end.
Could in principle take NTuples of arbitrary length by removing "3" from the type above.
"""
function TupleSmash(n1::NTuple{3,Int64}, n2::NTuple{3,Int64}, n3::NTuple{3,Int64}, level::Int64)
  return tuple(n1...,n2...,n3...,level)
end


"""
`EnumUp(nsum::Int64; fixed::Bool = false)`

Takes a dictionary argument and returns a list of all states possible from that one.
Lists all UP TO states with that many elements.  If fixed = true, returns elements which have 
elements summing to nsum ONLY.  

Test it:
find( x->isequal(x,(3,0,0)), outp) # replace with any tuple for (3,0,0)
find( x->isequal(x,(4,1,0)), outp)
OR:
length(outp)==length(unique(outp)) # works up to 25
for i = 0:25
  outp = EnumUp(i);
  println(length(outp)==length(unique(outp)))
end
OR: 
EnumUp(3; fixed = true) # will return 10 elements, instead of 20.  
"""
function EnumUp(nsum::Int64; fixed::Bool = false)
  termflag::Bool = true
  outp::Array{NTuple{3, Int64},1} = Array{NTuple{3,Int64},1}()
  push!(outp, (0,0,0))
  strt = 1
  if nsum > 0
    while termflag
      lng = length(outp)
      if lng > 1
        for el in strt+1:lng # won't doublecount, but won't work on first round.
          for nt in EnumerLevel(outp[el])
            push!(outp, nt)
          end
        end
      else # this should be the length 1 case only.
        for nt in EnumerLevel(outp[1])
          push!(outp, nt)
        end
        strt = length(outp)
      end
      if sum(outp[end]) == nsum
        termflag = false
      end
      strt = lng # reassign.
    end
  end
  if !fixed
    return outp
  else 
    return filter(x->sum(x)==nsum, outp) # this should only return elements with exactly the sum nsum.
  end 
end



"""
`HospitalDemand(pats::ProjectModule.patientcollection, )`
Takes a patientcollection and a whole state.
Makes a dict of every hospital's patients.
#TODO - does not take level into consideration.  FIX

Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
patients = CreateZips(ProjectModule.alldists, Texas);
FillPPatients(patients , ProjectModule.pinsured); # this line probably not necessary.
hzips, privrange, medrange,  unfound = HospitalDemand(patients);

"""
function HospitalDemand(pats::ProjectModule.patientcollection; 
                        mat1::Array{Float64,2} = ProjectModule.pcount)
  inter::Dict{Int64, Array{Int64,1}} = Dict{Int64, Array{Int64,1}}()
  # These must have the fid keys added to them.  
  outp_p::Dict{Int64, ProjectModule.patientrange} = Dict{Int64, ProjectModule.patientrange}()
  outp_m::Dict{Int64, ProjectModule.patientrange} = Dict{Int64, ProjectModule.patientrange}()
  # creates the objects to be mapped back
  for zp in keys(pats.zips) # these are zips
    for fid in keys(pats.zips[zp].facilities) # these are hospital fids
      if haskey(inter, fid)
        push!(inter[fid], zp)
      else 
        inter[fid] = Array{Int64,1}()
        push!(inter[fid], zp)
      end 
    end 
  end 
  # do some pre-processing on the pcount matrix - can write each zip as a dict
  d_p_l::Dict{Int64, patientcount} = Dict{Int64, patientcount}() # private, low 
  d_p_h::Dict{Int64, patientcount} = Dict{Int64, patientcount}() # private, high
  d_m_l::Dict{Int64, patientcount} = Dict{Int64, patientcount}() # medicaid, low 
  d_m_h::Dict{Int64, patientcount} = Dict{Int64, patientcount}() # medicaid, high
  for r in 1:7:size(mat1, 1) # skips by seven over DRGs
    d_p_l[mat1[r,1]] = patientcount(mat1[r,41],mat1[r+1,41],mat1[r+2,41],mat1[r+3,41],mat1[r+4,41],mat1[r+5,41],mat1[r+6,41]) # contains zipmin_p
    d_p_h[mat1[r,1]] = patientcount(mat1[r,40],mat1[r+1,40],mat1[r+2,40],mat1[r+3,40],mat1[r+4,40],mat1[r+5,40],mat1[r+6,40]) # contains zipmax_p
    d_m_l[mat1[r,1]] = patientcount(mat1[r,45],mat1[r+1,45],mat1[r+2,45],mat1[r+3,45],mat1[r+4,45],mat1[r+5,45],mat1[r+6,45]) # contains zipmin_m
    d_m_h[mat1[r,1]] = patientcount(mat1[r,44],mat1[r+1,44],mat1[r+2,44],mat1[r+3,44],mat1[r+4,44],mat1[r+5,44],mat1[r+6,44]) # contains zipmax_m
  end 

  # Now create the patientranges from these sets.  
  unfound = Set{Int64}()
  for k1 in keys(inter) # this is a Dict{Fid, Array{zip}}, so these are fids.
    if !haskey(outp_p, k1)
      outp_p[k1] = patientrange(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    end 
    if !haskey(outp_m, k1)
      outp_m[k1] = patientrange(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    end 
    for el in inter[k1] # this is an array of zip codes 
      if (haskey(d_p_l, el))&(haskey(d_p_h, el))&(haskey(d_m_l, el))&(haskey(d_m_h, el))
          # Privately insured.
        outp_p[k1].l385 += d_p_l[el].count385
        outp_p[k1].u385 += d_p_h[el].count385
        outp_p[k1].l386 += d_p_l[el].count386
        outp_p[k1].u386 += d_p_h[el].count386
        outp_p[k1].l387 += d_p_l[el].count387
        outp_p[k1].u387 += d_p_h[el].count387
        outp_p[k1].l388 += d_p_l[el].count388
        outp_p[k1].u388 += d_p_h[el].count388
        outp_p[k1].l389 += d_p_l[el].count389
        outp_p[k1].u389 += d_p_h[el].count389
        outp_p[k1].l390 += d_p_l[el].count390
        outp_p[k1].u390 += d_p_h[el].count390
        outp_p[k1].l391 += d_p_l[el].count391
        outp_p[k1].u391 += d_p_h[el].count391
          # Medicaid Patients 
        outp_m[k1].l385 += d_m_l[el].count385
        outp_m[k1].u385 += d_m_h[el].count385
        outp_m[k1].l386 += d_m_l[el].count386
        outp_m[k1].u386 += d_m_h[el].count386
        outp_m[k1].l387 += d_m_l[el].count387
        outp_m[k1].u387 += d_m_h[el].count387
        outp_m[k1].l388 += d_m_l[el].count388
        outp_m[k1].u388 += d_m_h[el].count388
        outp_m[k1].l389 += d_m_l[el].count389
        outp_m[k1].u389 += d_m_h[el].count389
        outp_m[k1].l390 += d_m_l[el].count390
        outp_m[k1].u390 += d_m_h[el].count390
        outp_m[k1].l391 += d_m_l[el].count391
        outp_m[k1].u391 += d_m_h[el].count391
      else
        push!(unfound,el) 
      end 
    end 
  end 
  return inter, outp_p, outp_m, unfound
end 






### Above this line... Exact Value development ###



"""
`ValApprox(D::DynState)`
This computes the dynamic simulation across all of the facilities in all of the markets.

TODO - where is the shock added and what is the variance of the shock?  It should be the normalizing constant.

To start:
dyn = CounterObjects();
V = allvisits(Dict{Int64, vrecord}());
ValApprox(dyn, V, 100_000 ; chunk = [2]) # just doing one hospital.
22.001426 seconds (281.50 M allocations: 5.099 GiB, 6.99% gc time)
"""
function ValApprox(D::DynState, V::allvisits, itlim::Int64; chunk::Array{Int64,1} = collect(1:size(D.all,1)), debug::Bool = false)
  iterations::Int64 = 0
  converged::Bool = false
  a::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  b::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  # NB - the next line is not necessary while we aren't doing oblivious yet.
  #steadylevs = AllAgg(D, chunk)
  for el in chunk                                                          # creates a dictionary of visited records.
    V.all[D.all[el].fid] = vrecord( Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, 1}(), 1)
  end
  while (iterations<itlim)&&(!converged)
    for el in D.all[chunk]
      if !el.converged                                                     # only keep simulating with the ones which haven't converged
        DSimNew(el.mk, el.fid, a, b)
        GetProb(el)                                                        # this chooses the action by the other firms
        # TODO: add a market size check here to do something more like oblivious in large markets.
        if !haskey(el.visited, KeyCreate(el.cns, el.level))
          println("adding new counter")
          el.visited[KeyCreate(el.cns, el.level)]=nlrec(MD(ChoicesAvailable(el), StartingVals(el, a, b)), vcat(ChoicesAvailable(el),transpose(PolicyUpdate(StartingVals(el, a, b)))), Dict(k => 0 for k in ChoicesAvailable(el)) )
        end
        Action = ChooseAction(el)                                              # Takes an action and returns it.
        ComputeR(el, a, b, Action, iterations; debug = debug)                  # Computes the return to the action
        level::Int64 = LevelFunction(el, Action)                                       # Level may change with action, but for next period.
        if iterations <= 1_000_000
          push!(V.all[el.fid].visited, RTuple(el, Action))                     # Record the first million state-action pairs in a vector
        elseif iterations >1_000_000
          V.all[el.fid].visited[iterations%1_000_000] = RTuple(el, Action)     # Once this is a million entries long, start overwriting to keep track of only 1_000_000
        end
        el.previous = el.level                                              # Reassign current level to previous.
        el.level = level                                                    # Reassign current level, if it has changed or not.
        ExCheck(el)                                                         # Checks for exit
        FixNN(el)                                                           # Fixes the firms neighbors.
        iterations += 1                                                     # Update iteration count - TODO: delete after debugging.
        V.all[el.fid].totalcnt += 1                                         # Update the iteration count within the visit records.
        PatientZero(a,b) # resets both patientcounts to zero.
      end
      # TODO - uncomment convergence test when that is debugged.
      if iterations%1_000_000 == 0                                        # Check for convergence every million iterations
        CheckConvergence(el) # FIXME - this is the wrong signature for this function.
      end
    end
  end
  converged = Halt(D, chunk)                                                # Check to see if all firms in "chunk" have converged, then halt if they have.
end

"""
`KeytoTuple(x::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64})`
The visit record is a giant tuple - this cuts that into (neighbors, level) and (action).
The former is a tuple, the latter a scalar.
"""
function KeytoTuple(x::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64})
  return x[1:end-1], x[end]
end




"""
`CheckConvergence(h::simh, V::Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},1}; draws::Int64 = 100, demands::Int64 = 10, disc::Float64 = 0.95, debug::Bool = true)`
Check the convergence criterion in Collard-Wexler or Pakes McGuire.
This can be checked every million iterations.  When that happens,
reset the counters for every state in "visited".
"Visited" accessed by V.all[fid].visited
- Convergence should be assessed at *visited* states, something like unique(V), not all of them.
"""
function CheckConvergence(h::simh, V::Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},1}; draws::Int64 = 1000, disc::Float64 = 0.95, debug::Bool = true)
  outp::Array{Float64,1} = Array{Float64, 1}()
  pairs::Array{Tuple{Float64,Float64},1} = Array{Tuple{Float64, Float64},1}()
  totvisits::Int64 = 0
  a::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  b::ProjectModule.patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  states::Dict{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, Tuple{Float64,Float64}} = Dict{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, Tuple{Float64,Float64}}()
  itercount::Int64 = 0
  for k in unique(V)                                                                                              # only check this set of values visited in the last million iterations.
   k1::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64}, k2::Int64 = KeytoTuple(k)
   approxim::Float64 = 0.0
   for d = 1:draws
      origlevel::Int64 = h.level                                                                          # keep track of the level inside of the loop so that it can be reset.
      nextact::Int64 = convert(Int64, sample(h.visited[k1].psi[1,:], WeightVec(h.visited[k1].psi[2,:])))  # Take an action.  NB: LevelFunction takes Int64 argument in second place.
      h.level = LevelFunction(h, nextact)                                                                 # this level must be updated so that the profit computation is correct.
      DSimNew(h.mk, h.fid, a, b)
      # currdem::Tuple{ProjectModule.patientcount,ProjectModule.patientcount} = SimpleDemand(dems, h.level) # draw from the limited demand set.
      currpi::Float64 = SinglePay(h, a, b, nextact)                                              # Current period return, excluding continuation value.
      contval::Float64 = 0.0
      if haskey(h.visited, KeyCreate(h.cns, h.level))                                                     # check neighbors/level pair
        #FIXME - note here: ContError should not be present upon exit.
        contval = disc*WProb(h.visited[KeyCreate(h.cns, h.level)])
        if nextact != 11
          contval += disc*(ContError(h.visited[KeyCreate(h.cns, h.level)]))
        end
      else
        #FIXME - what is happening here now on "not available"  contval stays 0 - but that's going to make the error larger.
        # But these will still be 0 since that's what StartingVals is giving.
        println("Added entry")
        # Choices available?  Is that going to get the level right?  What about act?
        h.visited[KeyCreate(h.cns, h.level)]=nlrec(MD(ChoicesAvailable(h), StartingVals(h, currdem[1], currdem[2])), vcat(ChoicesAvailable(h),transpose(PolicyUpdate(StartingVals(h, currdem[1], currdem[2])))), Dict(k => 0 for k in ChoicesAvailable(h)) )
        if nextat!=11
          contval += disc*()
        end
      end
      h.level = origlevel                                                                                 # reset the level to the original value.
      approxim += (currpi+contval)                                                                        # this needs to be weighted by the right count
      PatientZero(a,b) # resets both patientcounts to zero.
   end
   push!(outp, (approxim/draws - h.visited[k1].aw[k2])^2)                                                 # TODO - replace this with a sum when confidence is reached in the outcome..
   push!(pairs, (approxim/draws, h.visited[k1].aw[k2]))
   states[KeyCreate(h.cns, h.level)] = (approxim/draws, h.visited[k1].aw[k2])
  end
  # if !debug
  #   for k1 in keys(h.visited)
  #     for k2 in keys(h.visited[k1].counter)
  #       itercount += h.visited[k1].counter[k2]                                                                 # how many iterations were made in total?
  #       h.visited[k1].counter[k2] = 1                                                                          # Reset all counter keys to 1 to track which states visited in last million iterations.
  #     end
  #   end
  # end
  return outp, itercount, states, pairs                                                                             # TODO: eventually divide former by latter. And drop states.
end




"""
`LogitCheck(h::simh)`
Check the estimated probabilities from the record of visits against the probabilities
from the original logit model.
"""
function LogitCheck(h::simh)
  for k in keys(h.visited)
    levl = (-1, -1)
    choices = [0]
    if k[end] == 1
      levl = (0,0)
    elseif k[end] == 2
      levl = (1,0)
    elseif k[end] == 3
      levl = (0,1)
    end
    if k[end] == 1
      choices = [10 2 1 11]
    elseif k[end] == 2
      choices = [5 10 6 11]
    elseif k[end] == 3
      choices =  [4 3 10 11]
    end
    levels = MktSize(neighbors(k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9]))
    prs = logitest(levl, levels[1], levels[2], levels[3], [h.cns.level105; h.cns.level205; h.cns.level305; h.cns.level1515; h.cns.level2515; h.cns.level3515; h.cns.level11525; h.cns.level21525; h.cns.level31525 ] )
    println("*****************")
    println("Neighbors: ", k[1:end-1])
    println("Level: ", k[end])
    println("Equilibrium Probabilities: ", choices,"  ", prs)
    println("Estimated Probabilities: ", h.visited[k].psi)
  end
end



"""
`AllDems(d::DynState, ch::Array{Int64,1})`
Generate a dictionary of tuples of demand states for every facility in the array `ch`.
"""
function AllDems(d::DynState, ch::Array{Int64, 1}; repcount::Int64 = 10)
  outp::Dict{Int64, Tuple{Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}}}} = Dict{Int64, Tuple{Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}}}}()
  for el in ch
    outp[d.all[el].fid] = ThreeDemands(d.all[el], repcount)
  end
  return outp
end



"""
`Halt(D::DynState)`
For each firm in the simulation, check whether it has converged or not.  Return false
when one firm has been found which hasn't.
"""
function Halt(D::DynState, chunk::Array{Int64,1})
  b::Bool = true
  for el in D.all[chunk]
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

TODO - this does not actually update the value of the int action.

"""
function ChooseAction(h::simh)
  act::Int64 = sample(h.visited[KeyCreate(h.cns, h.level)].psi[1,:], WeightVec(h.visited[KeyCreate(h.cns, h.level)].psi[2,:]))
  if act == 11
    h.exit = true
  end
  return act
end


"""
`ExCheck(h::simh)`
Check if the firm exited, then restart.
When the firm exits we also reset all of the neighbors.
"""
function ExCheck(h::simh)
  if (h.exit)|(h.level== -999)
    h.exit = false
    h.level = h.actual
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
  for el in keys(h.visited)
    println("**************************")
    println(el)
    for act in keys(h.visited[el].aw)
      println("Action: ", act, " Value: ", h.visited[el].aw[act] ,"  Probability: ", h.visited[el].psi[2,findin(h.visited[el].psi[1,:], act)][1], "  Frequency: ", h.visited[el].counter[act],)
    end
  end
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
An easier way to work with demand.
This will:
- Compute demand at some level.
- Update to another.
- Compute at that level.
- Compute at the third level.
- Return a tuple of demands for medicaid and private patients.
- Random draws can be made from those demands easily.
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
  h.level = h.actual #reassign the levels
  UpdateDUtil(h)
  FixNS(h)
  FixMainNs(h)
  return outp1, outp2, outp3
end


"""
`SimpleDemand(ds::Tuple{Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}}}, lev::Int64)`
Takes a gigantic tuple of tuples of arrays of patientcounts and returns a random element
depending on the level.
The input to this function is the output of the function `ThreeDemands(h::simh, totl::Int64)`
"""
function SimpleDemand(ds::Tuple{Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}},Tuple{Array{ProjectModule.patientcount,1},Array{ProjectModule.patientcount,1}}}, lev::Int64)
  if lev == 1
    return EasyDemand(ds[1])
  elseif lev == 2
    return EasyDemand(ds[2])
  elseif lev == 3
    return EasyDemand(ds[2])
  else #exited - return 0 demand.
    return emp::Tuple{ProjectModule.patientcount,ProjectModule.patientcount} = (ProjectModule.patientcount(0,0,0,0,0,0,0),ProjectModule.patientcount(0,0,0,0,0,0,0))
  end
end



"""
`AvgAggState(h::simh)`
To make the computation faster, compute the average state in the rest of the market.  Draw actions for
all of the neighbors, compute the aggregate state, then take the mean over such states.
Note that this will also reset the market to its original state before drawing actions for the neighbors
in every round.
This function is called in `AllAgg.`
"""
function AvgAggState(h::simh, drws::Int64)
  tem::Array{neighbors, 1} = Array{neighbors, 1}()
  interim::neighbors = neighbors(0,0,0,0,0,0,0,0,0)
  for i = 1:drws
    FixNS(h)
    GetProb(h)
    interim = sum(interim, h.cns)
  end
  h.level = h.actual #reset firm level to the actual
  FixNS(h)
  UpdateDUtil(h)
  FixMainNs(h)
  return neighbors( round(Int64, interim.level105/drws), round(Int64,interim.level205/drws ), round(Int64, interim.level305/drws), round(Int64, interim.level1515/drws), round(Int64, interim.level2515/drws), round(Int64, interim.level3515/drws), round(Int64,interim.level11525/drws ), round(Int64, interim.level21525/drws), round(Int64, interim.level31525/drws) )
end


"""
`AllAgg(d::DynState, ch::Array{Int64, 1})`
Computes three aggregate market states for each of the firms which are in "ch"
The return type is a dictionary of Tuple{Int64, Int64} keys mapping to a neighbors type.
"""
function AllAgg(d::DynState, ch::Array{Int64,1})
  outp::Dict{Tuple{Int64, Int64}, ProjectModule.neighbors} = Dict{Tuple{Int64, Int64}, ProjectModule.neighbors}()
  for el in d.all[ch]
    for levs = [1,2,3]
      FixNS(el)
      el.level = levs
      UpdateDUtil(el)
      outp[(el.fid, levs)] = AvgAggState(el, 20)
    end
    el.level = el.actual #reset firm level to the actual
    FixNS(el)
    UpdateDUtil(el)
    FixMainNs(el)
  end
  return outp
end

"""
`FixNS(h::simh)`
Take all of the neighbors and reset their states to their true levels if those levels aren't what they
should be
"""
function FixNS(h::simh)
  for nb in h.ns
    if nb.level != nb.truelevel # if the firm exited, reset it.
      nb.level = nb.truelevel
    end
  end
end




"""
`NeighborsCheck(d::DynState)`
How many neighbors do the firms have?
"""
function NeighborsCheck(d::DynState)
  for el in d.all
    println(el.fid, "  ", sum(el.cns))
  end
end


"""
`VisitPrint(h::simh)`
Print the states visited, the counts and the vals approximated
"""
function VisitPrint(h::simh)
  for el in keys(h.visited)
    println("State: ", el)
    println("Values: ", h.visited[el].aw)
    for k in keys(h.visited[el].counter)
      println("Action ", k, "  Counter ", h.visited[el].counter[k])
    end
  end
end

"""
`ProbChange(h::simh)`
Watch how the probabilities change
"""
function ProbChange(h::simh)
  println(h.visited[KeyCreate(h.cns, h.level)].psi)
end

"""
`findscale(d::DynState)`
Run a demand and profit computation a bunch of times to find the max to use as a scaling factor.
Also return the min just to get an idea of where the bottom is.
"""
function FindScale(d::DynState)
  maxprof::Float64 = 0.0
  minprof::Float64 = 1_000_000
  for el in d.all
    for i = 1:10
      a, b = DSim(el.mk, el.fid)
      if SinglePay(el, a, b) > maxprof
        maxprof = SinglePay(el, a, b)
      elseif SinglePay(el,a,b) < minprof
        minprof = SinglePay(el, a, b)
      end
    end
  end
  return maxprof, minprof
end

"""
`MktShare(cp::cpats)`
Compute the market share
"""

function MktShare(cp::cpats, fid::Int64)



end









#=


=#
