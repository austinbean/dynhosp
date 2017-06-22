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
  return DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);
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

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);


# another test - make a hospital first:
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);
dm, dp = PRanges(ProjectModule.pcount)
h1 = simh(4530190, 30.289991, -97.726196, 3, 3, 3, 100, neighbors(0,0,0,0,0,0,0,0,0), Array{Int64,1}(), Dict{Tuple{Int64}, ProjectModule.nlrec}(), Array{ProjectModule.shortrec,1}(), DynPatients(patients, 4530190, dm, dp), false, false, false)

Now use this to continue debugging.  Next is... 


Note that the function takes TWO EntireState arguments.  This is super dumb, but
only one of them (containing hospital types) has the bed counts.
"""
function DynStateCreate( Tex::EntireState, Tex2::EntireState, p::patientcollection, zipdrg::Array{Float64,2} )
  outp = DynState(Array{simh,1}())
  d_m, d_p = PRanges(zipdrg)
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
                     DynPatients(p, Tex.mkts[k1].collection[hk].fid, d_m, d_p), # this is the cmkt (?) should create the patient collection as a subelement of the hospital record.
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
    # Can be replaced with MakeNL 
    el.visited[KeyCreate(el.cns, el.level)] = MakeNL(el, ChoicesAvailable(el), ProjectModule.patientcount(5,6,4,13,8,41,248), ProjectModule.patientcount(5,6,4,13,8,41,248))
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
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);
d_m, d_p = PRanges(ProjectModule.pcount)

DynPatients(patients, 4530190, d_m, d_p);

"""
function DynPatients(p::patientcollection, f::Int64, d_m::Dict{Tuple{Int64,Int64},Tuple{Int64,Int64}} , d_p::Dict{Tuple{Int64,Int64},Tuple{Int64,Int64}} )
  outp::cmkt = cmkt(f, Array{cpats,1}())
  zpc = PatientFind(p, f) # finds the zip codes
  for el in zpc
    if haskey(d_m, (el,385))&haskey(d_p, (el,385)) # FIXME - think about those not found.
      push!(outp.m, cpats(el,
                    p.zips[el].lat,
                    p.zips[el].long,
                    DetUtils(p.zips[el]; switch = false),
                    DetUtils(p.zips[el]; switch = true),
                    vcat(transpose(DetUtils(p.zips[el]; switch = false)[1,:]), transpose(CounterWTP(DetUtils(p.zips[el]; switch = false)[2,:]))), #NB: bottom row only.
                    Array{shortrec,1}(),  
                    patientrange(d_m[(el,385)][1],d_m[(el,385)][2],d_m[(el,386)][1],d_m[(el,386)][2],d_m[(el,387)][1],d_m[(el,387)][2],d_m[(el,388)][1],d_m[(el,388)][2],d_m[(el,389)][1],d_m[(el,389)][2],d_m[(el,390)][1],d_m[(el,390)][2],d_m[(el,391)][1],d_m[(el,391)][2]), 
                    patientrange(d_p[(el,385)][1],d_p[(el,385)][2],d_p[(el,386)][1],d_p[(el,386)][2],d_p[(el,387)][1],d_p[(el,387)][2],d_p[(el,388)][1],d_p[(el,388)][2],d_p[(el,389)][1],d_p[(el,389)][2],d_p[(el,390)][1],d_p[(el,390)][2],d_p[(el,391)][1],d_p[(el,391)][2]) ) )  
    else 
      # println(el) # lists zips which are not found.  
    end 
  end
  return outp
end


"""
`PRanges(p1::Array{Float64,2})`

This will return two dicts.  Each dict is d[(zipcode, drg)] = (max, min)
where this both the key and the return are tuples.

d_m, d_p = PRanges(ProjectModule.pcount)

"""
function PRanges(p1::Array{Float64, 2})
  maxes_m::Int64 = 44 # Medicaid patients 
  mins_m::Int64 = 45
  maxes_p::Int64 = 40 # Private patients. 
  mins_p::Int64 = 41
  zips::Int64 = 1
  drgs::Int64 = 2
  outp_m::Dict{Tuple{Int64, Int64}, Tuple{Int64,Int64}} = Dict{Tuple{Int64,Int64}, Tuple{Int64,Int64}}()
  outp_p::Dict{Tuple{Int64, Int64}, Tuple{Int64,Int64}} = Dict{Tuple{Int64,Int64}, Tuple{Int64,Int64}}()
  for i = 1:size(p1,1) #rows 
    for j = 1:size(p1,2) #columns
      outp_m[(convert(Int64,p1[i,zips]),convert(Int64,p1[i,drgs]))] = (convert(Int64,p1[i,mins_m]),convert(Int64,p1[i,maxes_m]))
      outp_p[(convert(Int64,p1[i,zips]),convert(Int64,p1[i,drgs]))] = (convert(Int64,p1[i,mins_p]),convert(Int64,p1[i,maxes_p]))
    end 
  end 
  return outp_m, outp_p
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
  # FIXME - maybe start this at 1.0 to deal with infinity problem?
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


### Testing ###
temparr = zeros(2, 12)
WTPNew(dyn.all[2].mk.m[1].putils, temparr)
dyn.all[2].mk.m[1].putils[2,:] = zeros(1,6)
ArrayZero(temparr)

# prior version (not correcting for occasional NaN's) 0.000001 / 0 allocations.
# new version: 0.000001 seconds (3 allocations: 208 bytes)
for i = 1:10
  new_r = rand(2,10)
  @time WTPNew(new_r, temparr)
  ArrayZero(temparr)
end 


"""
function WTPNew(c::Array{Float64,2}, arr::Array{Float64,2})
  # doing int_sum = 1.0 captures the fact that there is one outside option, with util 0, so int_sum must
  # be at least one.  This solves the division by zero problem.  plus all elements which are zero can be ignored.
  # the first zero in the top row receives 1/int_sum in the second row.  This gets the prob of the OO.   
  int_sum::Float64 = 1.0 
  for el in 1:size(c,2)
    if c[1,el]!=0.0 
      arr[1,el] = c[1,el]
      arr[2,el] = exp(c[2,el])
      int_sum += arr[2,el]
    end
  end
  for i=1:size(c,2)
    arr[2,i]/=int_sum   
  end
  arr[2,findfirst(arr[2,:], 0)] = 1/int_sum # this is a little expensive, but not too bad.
end


"""
`DemComp(inparr::Array{Float64,2}, fid::Int64, c::cpats, p_or_m::Bool)`
Takes a fid.  Finds the index in inparr.  Takes the share from inparr.
Multiplies all of patient values in cpat by that number.
Returns the value.


TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);;

p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
ptest = patientcount(1,2,3,4,5,6,7)
temparr = zeros(2, 12) # 12 is max # of facilities. 


DemComp(dyn.all[2].mk.m[1].putils, temparr, p1, dyn.all[2].fid,   ptest)

function TestDemComp(n::Int64)
  p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

  temparr = zeros(2, 12) # 12 is max # of facilities. 

  for i = 1:n
    @time DemComp(dyn.all[2].mk.m[1].putils, temparr, p1, dyn.all[2].fid, dyn.all[2].mk.m[1], true, false )
    PatientZero(p1, p2)
  end 

end 

TestDemComp(5)
"""
function DemComp(inparr::Array{Float64,2}, temparr::Array{Float64,2}, pp::patientcount,fid::Int64, c::patientcount)  
  # NB: inparr is a sub-field of cpats.  inparr is either c.putils or c.mutils.
  index::Int64 = 0
  counter::Int64 = 0 # Records the index  ie.  p385, p386, etc
  for i = 1:size(inparr,2)
    if inparr[1,i] == fid
      index = i #reassign - XXX - double check that this reassignment works.  I think it does.  
    end
  end
  WTPNew(inparr, temparr) # updates temparr
  if index!=0 # don't look for a facility that isn't there.
    for j in c
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
dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);;



p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
DSimNew(dyn.all[2].mk, dyn.all[2].fid, p1, p2)

function TestDSimNew(n::Int64)
  p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  for i = 1:n 
    @time DSimNew(dyn.all[2].mk, dyn.all[2].fid, p1, p2)
    #println(p1, "  ", p2)
    PatientZero(p1, p2)
  end
end 

TestDSimNew(5)

"""
function DSimNew(c::cmkt, f::Int64, pcount::patientcount, mcount::patientcount; maxh::Int64 = 12)
  temparr::Array{Float64,2} = zeros(2, maxh)
  for el in c.m
    ArrayZero(temparr)
    DemComp(el.putils, temparr, pcount, f, DrawAll(el.pcounts))  # pcount is an empty, pre-allocated patientcount into which results are written.
    DemComp(el.mutils, temparr, mcount, f, DrawAll(el.mcounts)) 
  end
end


"""
`UpdateDUtil(h::simh)`
When something in the market changes, the utility must be updated in all zip codes
for which it can be chosen.  This calls the function HUtil to do the actual update.

##  To test ## 
  # one test: insert these lines after h.tbu - will check that coordinate found is correct.  
  # This test passes.
      if !isapprox(el.putils[1,findin(el.putils[1,:], h.fid)][1], h.fid) 
        println("private doesn't match: ", el.putils[1,findin(el.putils[1,:], h.fid)][1])
        println("fid: ", h.fid )
      end 
      if !isapprox(el.mutils[1,findin(el.mutils[1,:], h.fid)][1], h.fid) 
        println("medicaid doesn't match: ", el.mutils[1,findin(el.mutils[1,:], h.fid)][1])
        println("fid: ", h.fid )
      end 


"""
function UpdateDUtil(h::simh)
  for el in h.mk.m
    for sr in el.facs
      if sr.tbu # if this is true, the hospital needs updating
        # this is a collection of shortrecs within a zip.  
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

# FIXME - this needs to include the possibility of the shock leading to a different choice.  
"""
function GetProb(s::simh)
  for el in s.ns
    if el.level != -999
      # FIXME - here is the place to fix the variance problem, I think.  
      action::Int64 = sample(ChoicesAvailable(el), el.chprobs )                            # Take the action
      di::Float64 = distance(el.lat, el.long, s.lat, s.long)
      if action == 1
        el.level = 3
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs = Weights(vec(logitest((0,1), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 2
        el.level = 2
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =Weights(vec(logitest((1,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 3
        el.level = 2
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =Weights(vec(logitest((1,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 4
        el.level = 1
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =Weights(vec(logitest((0,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 5
        el.level = 1
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =Weights(vec(logitest((0,0), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
        el.tbu = true
      elseif action == 6
        el.level = 3
        el.choices = ChoicesAvailable(el)
        levels = MktSize(el.ns)
        el.chprobs =Weights(vec(logitest((0,1), levels[1], levels[2], levels[3], [el.ns.level105; el.ns.level205; el.ns.level305; el.ns.level1515; el.ns.level2515; el.ns.level3515; el.ns.level11525; el.ns.level21525; el.ns.level31525 ])))
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
`StartingVals(h::simh)`
This will generate a set of starting values for the hospitals in the market.
Idea - just take the PDV of one period's profit.  But it needs to be done
over the whole set of possible actions.  Right now this assumes that *all* of the
actions have the same value, except for exit.  That can be improved.

This is pretty dumb.   
"""
function StartingVals(h::simh,
                      ppats::patientcount,
                      mpats::patientcount;
                      disc::Float64 = 0.95)
       #return [0.0,0.0,0.0,0.0]
  return vcat(repmat([min(SinglePay(h, ppats, mpats, 10)/(1-disc), 100.0)],3), [min(SinglePay(h, ppats, mpats, 10)/((1-disc)*1000), 100.0)])
end




"""
`ProbUpdate(aw::Dict{Int64,Float64})`
This should update the probabilities.


for i = 1:10
  @time ProbUpdate(aw, arr)
end 

# FIXME - is there a good reason to do this in two functions?  Maybe used elsewhere? 
ProbUpdate(dyn.all[1].visited[(0,0,0,0,0,0,1,0,0,1)].aw, dyn.all[1].visited[(0,0,0,0,0,0,1,0,0,1)].psi)

"""
function ProbUpdate(aw::Dict{Int64,Float64}, outp::Array{Float64,2})
  # outp::Array{Float64,2} = Array{Float64,2}(2,4)
  for (ind,el) in enumerate(keys(aw))
    outp[1,ind] = el 
    outp[2,ind] = aw[el]
  end
  PolicyUpdate(outp)  
  return outp   
end


"""
`PolicyUpdate(neww::Array{Float64,1}; ep::Float64 = 0.0001)`
Takes an array of the new W values and maps back to the simh action probabilities.
NB: When probs are 0, continuation value of error is problematic, since this is
given by the log, so the value is constrained to be this small positive value.
The problem is basically underflow: when returns are really high to staying in and
much smaller to getting out, the estimated prob is zero.

Test this:

for i = 1:10
  ab = rand(2,4)
  @time PolicyUpdate(ab)
end 

0.000001 seconds (5 allocations: 368 bytes)
"""
function PolicyUpdate(neww::Array{Float64,2}; ep::Float64 = 0.000000001)
  mx::Float64 = maximum(neww[2,:])
  sm::Float64 = sum(exp.(neww[2,:] .- mx))
  for i = 1:size(neww, 2)
    neww[2,i] = max( (exp(neww[2,i]-mx))/sm, ep) 
  end 
end


"""
`PolicyUp2(newch::Array{Int64,2}, neww::Array{Float64,1}; ep::Float64 = 0.000000001)`
Just in the case where the record is initialized we need a different input category for this function.
This is pretty dumb, but it would slow down PolicyUpdate to fix it and this function isn't called often,
only when records need to be initialized in DynStateCreate.

Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);
dm, dp = PRanges(ProjectModule.pcount)
h1 = simh(4530190, 30.289991, -97.726196, 3, 3, 3, 100, neighbors(0,0,0,0,0,0,0,0,0), Array{Int64,1}(), Dict{Tuple{Int64}, ProjectModule.nlrec}(), Array{ProjectModule.shortrec,1}(), DynPatients(patients, 4530190, dm, dp), false, false, false)

ch1 = ChoicesAvailable(h1)
pr1 = StartingVals(h1, patientcount(1,1,1,1,1,1,1), patientcount(1,1,1,1,1,1,1,))

pr2 = [ 10.0, 10.0, 0.0, 0.0]
PolicyUp2(ch1, pr1)
PolicyUp2(ch1, pr2)
"""
function PolicyUp2(newch::Array{Int64,2}, neww::Array{Float64,1}; ep::Float64 = 0.000000001)
  outp::Array{Float64,2} = zeros(2, size(newch, 2))
  mx::Float64 = maximum(neww)
  sm::Float64 = sum(exp.(neww.-mx))
  for i = 1:size(newch,2)
    outp[1,i] = newch[i]
    outp[2,i] = max(exp(neww[i]-mx)/sm, ep)
  end 
  return outp 
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
  nmes::Array{Int64,1}=Array{Int64,1}()
  outp::Array{Float64,1}=Array{Float64,1}()
  for k in keys(d)
    push!(nmes, k)
    push!(outp, d[k])
  end
  return transpose(hcat(nmes,  PolicyUpdate(outp)))
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
  return (eulergamma - dot(log.(n.psi[2,:]), n.psi[2,:])) 
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

### Testing: ###
# Create hospital record:
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);
dm, dp = PRanges(ProjectModule.pcount)


h1 = simh(4530190, 30.289991, -97.726196, 3, 3, 3, 100, neighbors(0,0,0,0,0,0,0,0,0), Array{Int64,1}(), Dict{Tuple{Int64}, ProjectModule.nlrec}(), Array{ProjectModule.shortrec,1}(), DynPatients(patients, 4530190, dm, dp), false, false, false)
p1 = patientcount(100,100,100,100,100,100,100)
p2 = patientcount(100,100,100,100,100,100,100)


SinglePay(h1, p1, p2, 10) == 0.073280454324209 # returns true. 

"""
function SinglePay(s::simh,
                    mpats::ProjectModule.patientcount,
                    ppats::ProjectModule.patientcount,
                    action::Int64)
  # CONSTANTS:
    const scalefact::Float64 = 3.0e9
    const alf1::Float64 = 8336.17
    const alf2::Float64 = 36166.6
    const alf3::Float64 = 16309.47
    const gamma_1_385::Float64 = 20680.0 # ✓
    const gamma_2_385::Float64 = 42692.37 # ✓
    const gamma_3_385::Float64 = 20962.97 # ✓
    const gamma_1_386::Float64 = 81918.29 # X
    const gamma_2_386::Float64 = 74193.4 # X
    const gamma_3_386::Float64 = 99065.79 # X
    const gamma_1_387::Float64 = 30405.32 # X
    const gamma_2_387::Float64 = 49801.84 # X
    const gamma_3_387::Float64 = 22376.8 # X
    const gamma_1_388::Float64 = 10051.55 # ✓
    const gamma_2_388::Float64 = 19019.18 # X
    const gamma_3_388::Float64 = 33963.5 # X
    const gamma_1_389::Float64 = 29122.89 # X
    const gamma_2_389::Float64 = 14279.58 # X
    const gamma_3_389::Float64 = 20708.15 # X
    const gamma_1_390::Float64 = 22830.05 # X
    const gamma_2_390::Float64 = 6754.76 # X
    const gamma_3_390::Float64 = 3667.42 # ✓
    const gamma_1_391::Float64 = 9089.77 # X
    const gamma_2_391::Float64 = 8120.85 # X
    const gamma_3_391::Float64 = 1900.5 # ✓
    const level12::Float64 = 1.64669492e6
    const level13::Float64 = 5.0165876e6
    const level21::Float64 = -366430.33
    const level23::Float64 = 1.83969306e6
    const level31::Float64 = -90614.32
    const level32::Float64 = -157206.98
    const mcaid385::Float64 = 151380.0
    const mcaid386::Float64 = 48417.0
    const mcaid387::Float64 = 18845.0
    const mcaid388::Float64 = 7507.0
    const mcaid389::Float64 = 9424.0
    const mcaid390::Float64 = 4623.0
    const mcaid391::Float64 = 3664.0 # to DRG mean added 3094 - avg reimbursement for DRGs 370-375 under TX Medicaid (2012)
  # Compute WTP.
    outp::Float64 = 0.0
    levelc::Float64 = 0.0
    wtp::Float64 = FindWTP(s)
    if (s.level == 1)&(action!=11)
      outp = alf1*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_1_385*(ppats.count385+mpats.count385) - gamma_1_386*(ppats.count386+mpats.count386) - gamma_1_387*(ppats.count387+mpats.count387) - gamma_1_388*(mpats.count388+ppats.count388) - gamma_1_389*(mpats.count389+ppats.count389) - gamma_1_390*(ppats.count390+mpats.count390) - gamma_1_391*(ppats.count391+mpats.count391)
    elseif (s.level == 2)&(action!=11)
      outp = alf2*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_2_385*(ppats.count385+mpats.count385) - gamma_2_386*(ppats.count386+mpats.count386) - gamma_2_387*(ppats.count387+mpats.count387) - gamma_2_388*(mpats.count388+ppats.count388) - gamma_2_389*(mpats.count389+ppats.count389) - gamma_2_390*(ppats.count390+mpats.count390) - gamma_2_391*(ppats.count391+mpats.count391)
    elseif (s.level == 3)&(action!=11)
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
  # Current issue - return to changing much higher than return to staying, but that does not make sense...
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
    # FIXME - something to work on in ProbUpdate.  
    ProbUpdate(hosp.visited[k1].aw, hosp.visited[k1].psi) #WeightedProbUpdate(hosp.visited[k1].aw, hosp.visited[k1].psi, iterations)
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
`MakeNL(hosp::simh, chs::Array{Int64,2}, ppats::patientcount, mpats::patientcount)`
Returns nlrec - necessary?  
"""
function MakeNL(hosp::simh, chs::Array{Int64,2}, ppats::patientcount, mpats::patientcount)
  prbs::Array{Float64,2} = zeros(2,4)
  vals::Array{Float64,1} = StartingVals(hosp, ppats, mpats)
  d1::Dict{Int64,Int64} = Dict{Int64,Int64}()
  for i = 1:size(prbs,2)
    prbs[1,i] = chs[i]
    prbs[2,i] = vals[i]
    d1[chs[i]] = 0
  end 
  outp::nlrec = nlrec(MD(chs, vals), prbs, d1)
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
`StatePermute(d::DynState, states::Array{Tuple{Int64,Int64}})`
Takes a row of the state configurations.  
These are (competitor location in d.all,level) tuples.
Sets the d.all[location].level = level 
Sets d.all[location].tbu = true - should force update of DUtil and WTP.

"""
function StatePermute(d::DynState, states::Array{Tuple{Int64,Int64}})
  for (loc, lev) in states
    d.all[loc].level = lev 
    d.all[loc].tbu = true 
  end 
end 


"""
`StateBlock(n::Array{Tuple{Int64,Int64},2})`
Take what is given and append. 

ab = [ (10,1);
       (10,2);
       (10,3)]
returns: [ (20, 1),(10, 1);
           (20, 2),(10, 1);
           (20, 3),(10, 1);
           ... ]

StateBlock([(10,1); (10,2); (10,3); (10,999)], 20)
StateBlock([], 10)
"""
function StateBlock(n::Array, i::Int64) 
  len, width = size(n,1,2)
  if len > 1 # add these elements to whatever array exists.  
    outp::Array{Tuple{Int64,Int64},2} = Array{Tuple{Int64,Int64},2}(4*len, width+1)
    for el in 1:size(n,1) # rows of input 
        outp[4*(el-1)+1,1] = (i,1)
        outp[4*(el-1)+2,1] = (i,2)
        outp[4*(el-1)+3,1] = (i,3)
        outp[4*(el-1)+4,1] = (i,999)
      for (ix,tp) in enumerate(n[el,:]) # enumerate row elements.
        outp[4*(el-1)+1,ix+1] = tp 
        outp[4*(el-1)+2,ix+1] = tp 
        outp[4*(el-1)+3,ix+1] = tp
        outp[4*(el-1)+4,ix+1] = tp
      end   
    end 
    return outp 
  else # this is the base case.
    return [(i,1); (i,2); (i,3); (i,999)]
  end 
end 

"""
`MakeStateBlock(n::Array{Int64,1})`
Calls StateBlock repeatedly to create a list of all states. 
What this returns for a list of fids given in n is a collection of
possible organizations of the states of those firms.  E.g, this is called
on a list of neighbors n1, ..., nm and returns a list of tuples:
[ (n1, i), (n2, j), ..., (nm, k) ] where i, j, k ∈ {1,2,3,999} covering
different levels (999 is exit).

Not efficient due to reallocation of outp every time, but only 
needs to be called once, so not so bad.  

### Test ###

MakeStateBlock([10])
MakeStateBlock([10, 20])
MakeStateBlock([10, 20, 30])

"""
function MakeStateBlock(n::Array{Int64,1})
  outp::Array{Tuple{Int64,Int64}} = Array{Tuple{Int64,Int64}}()
  for i in n 
    outp = StateBlock(outp, i)
  end 
  return outp 
end 






"""
`ConvTest(d1::Dict{Int64,Bool})`
Iterates over the elements in totest and returns true when 
all are true, else false.


"""
function ConvTest(d1::Dict{Int64, Bool})
  outp::Bool = true 
  for k1 in keys(d1)
    outp = outp&d1[k1]
  end 
  return !outp # the while-loop in ExactVal continues when this returns true.  See tests above.  
end 





"""
`ContVal( futures::Dict{NTuple{9,Int64},Float64}, fid::Int64, stable::Dict{ Int64, Dict{NTuple{10, Int64},  Float64 } }, lev::Int64)`


Takes the states of other firms in the Dict futures, computes continuation vals using the Dict stable, 
then returns a float of the CV.  
futures should be the dict returned by TotalCombine - this is a set of probabilities of various outcomes.
stable is the dict of values - i.e., the continuation values of the firms
level is the level for the firm for which we are computing CV 
fid is the dyn.all[xxx].fid - the firm for which we are computing the CV.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);;


d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

d1[dyn.all[18].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[18].fid][StateKey(dyn.all[18], 1)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 2)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 3)] = 0.0
location2 = FindComps(dyn, dyn.all[18]) # locations are 19 and 152
push!(location2, 18)
# FIX THIS - above should be a dict. 
loc1 = Dict(672285 => 19, 670132 => 18, 373510 => 152 ) 
out_1 = Dict{Int64,NTuple{9,Int64}}()
StateRecord(loc1, dyn, out_1) 

d1[dyn.all[19].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[19].fid][TAddLevel(out_1[dyn.all[19].fid], 1)] = 0.0
d1[dyn.all[19].fid][TAddLevel(out_1[dyn.all[19].fid], 2)] = 0.0
d1[dyn.all[19].fid][TAddLevel(out_1[dyn.all[19].fid], 3)] = 0.0

d1[dyn.all[152].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[152].fid][TAddLevel(out_1[dyn.all[152].fid], 1)] = 0.0
d1[dyn.all[152].fid][TAddLevel(out_1[dyn.all[152].fid], 2)] = 1.0
d1[dyn.all[152].fid][TAddLevel(out_1[dyn.all[152].fid], 3)] = 2.0

cp2 = ContProbs(dyn.all[18].fid, out_1, d1)

compprobs = TotalCombine(dyn, 18, loc1, cp2)

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
      else 
        println("In ContVal, adding a state which should be there! ")
        stable[fid][TAddLevel(k1,lev)] = 0.5
      end 
    end 
  elseif lev == 2
    for k1 in keys(futures)
      if haskey(stable[fid],TAddLevel(k1,lev) )
        outp += futures[k1]*stable[fid][TAddLevel(k1,lev)]
      else 
        stable[fid][TAddLevel(k1,lev)] = 0.5
        println("In ContVal, adding a state which should be there! ")
      end     end 
  elseif lev == 3
    for k1 in keys(futures)
      if haskey(stable[fid],TAddLevel(k1,lev) )
        outp += futures[k1]*stable[fid][TAddLevel(k1,lev)]
      else 
        stable[fid][TAddLevel(k1,lev)] = 0.5
        println("In ContVal, adding a state which should be there! ")
      end 
    end
  else 
    #do nothing.  
  end 
  return outp 
end 





"""
`ContProbs(fid::Int64, state_recs::Dict{Int64,NTuple{9,Int64}}, stable_vals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64} })
A rewrite of the first ContProbs function. 
- (fid) Takes the FID of the main facility as an argument.
- (state_recs) State Recs of restricted states from the point of view of all neighbors.
- (nlocs) Locations of neighbors.  (necessary...)
- (stable_vals) The current equilibrium value estimations.
- (D) DynState.  (not necessary?)


### TESTING ###

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);; 
d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()

d1[dyn.all[11].fid] = Dict{NTuple{10,Int64}, Float64}() # this gives the argument stable_vals
d1[dyn.all[11].fid][StateKey(dyn.all[11], 1)] = 0.0
d1[dyn.all[11].fid][StateKey(dyn.all[11], 2)] = 0.0
d1[dyn.all[11].fid][StateKey(dyn.all[11], 3)] = 0.0

location2 = FindComps(dyn, dyn.all[11])

testcp2 = Dict{Int64,Int64}()

CompsDict(FindComps(dyn, dyn.all[11]), dyn, testcp2) # this is argument nlocs


StateRecord(testcp2, dyn, testcp22)

d1[dyn.all[location2[1]].fid] = Dict{NTuple{10,Int64}, Float64}(); d1[dyn.all[location2[2]].fid] = Dict{NTuple{10,Int64}, Float64}(); d1[dyn.all[location2[3]].fid] = Dict{NTuple{10,Int64}, Float64}(); d1[dyn.all[location2[4]].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[location2[1]].fid][NStateKey(testcp22[dyn.all[location2[1]].fid], 1)] = 0.0; d1[dyn.all[location2[1]].fid][NStateKey(testcp22[dyn.all[location2[1]].fid], 2)] = 0.0; d1[dyn.all[location2[1]].fid][NStateKey(testcp22[dyn.all[location2[1]].fid], 3)] = 0.0;
d1[dyn.all[location2[2]].fid][NStateKey(testcp22[dyn.all[location2[2]].fid], 1)] = 0.0; d1[dyn.all[location2[2]].fid][NStateKey(testcp22[dyn.all[location2[2]].fid], 2)] = 0.0; d1[dyn.all[location2[2]].fid][NStateKey(testcp22[dyn.all[location2[2]].fid], 3)] = 0.0;
d1[dyn.all[location2[3]].fid][NStateKey(testcp22[dyn.all[location2[3]].fid], 1)] = 0.0; d1[dyn.all[location2[3]].fid][NStateKey(testcp22[dyn.all[location2[3]].fid], 2)] = 0.0; d1[dyn.all[location2[3]].fid][NStateKey(testcp22[dyn.all[location2[3]].fid], 3)] = 0.0;
d1[dyn.all[location2[4]].fid][NStateKey(testcp22[dyn.all[location2[4]].fid], 1)] = 0.0; d1[dyn.all[location2[4]].fid][NStateKey(testcp22[dyn.all[location2[4]].fid], 2)] = 0.0; d1[dyn.all[location2[4]].fid][NStateKey(testcp22[dyn.all[location2[4]].fid], 3)] = 0.0;

ContProbs(dyn.all[11].fid, testcp22,  d1) # note that this takes the location of the main firm.  

# do it for entry 18 in dyn.  
d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d1[dyn.all[18].fid] = Dict{NTuple{10,Int64}, Float64}() # this gives the argument stable_vals
d1[dyn.all[18].fid][StateKey(dyn.all[18], 1)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 2)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 3)] = 0.0

location2 = FindComps(dyn, dyn.all[18])
testcp2 = Dict{Int64,Int64}()
CompsDict(FindComps(dyn, dyn.all[18]), dyn, testcp2) # this is argument nlocs
testcp22 = Dict{Int64, NTuple{9,Int64}}() # holds state_recs

StateRecord(testcp2, dyn, testcp22)
d1[dyn.all[location2[1]].fid] = Dict{NTuple{10,Int64}, Float64}(); d1[dyn.all[location2[2]].fid] = Dict{NTuple{10,Int64}, Float64}(); 
d1[dyn.all[location2[1]].fid][NStateKey(testcp22[dyn.all[location2[1]].fid], 1)] = 0.0; d1[dyn.all[location2[1]].fid][NStateKey(testcp22[dyn.all[location2[1]].fid], 2)] = 0.0; d1[dyn.all[location2[1]].fid][NStateKey(testcp22[dyn.all[location2[1]].fid], 3)] = 0.0;
d1[dyn.all[location2[2]].fid][NStateKey(testcp22[dyn.all[location2[2]].fid], 1)] = 0.0; d1[dyn.all[location2[2]].fid][NStateKey(testcp22[dyn.all[location2[2]].fid], 2)] = 0.0; d1[dyn.all[location2[2]].fid][NStateKey(testcp22[dyn.all[location2[2]].fid], 3)] = 0.0;


test1 = ContProbs(dyn.all[18].fid, testcp22, d1) 
# returns: Dict(672285 => [0.333333, 0.333333, 0.333333], 373510 => [0.333333, 0.333333, 0.333333])



dict1 = Dict('a'=> 1, 'b'=>2)
dict2 = Dict('a'=> 3, 'b'=>4)
keys(dict1) == keys(dict2)

"""
function ContProbs(fid::Int64,
             state_recs::Dict{Int64,NTuple{9,Int64}},
             stable_vals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64} })
  outp::Dict{Int64, Array{Float64,1}} = Dict{Int64, Array{Float64,1}}()
  for el in keys(stable_vals) # these should be all the fids.
    if el != fid  # we skip continuation probs for one firm.   
      outp[el] = exp.([stable_vals[el][TAddLevel(state_recs[el], 1)], stable_vals[el][TAddLevel(state_recs[el],2)], stable_vals[el][TAddLevel(state_recs[el], 3)]])
      outp[el] ./=(sum(outp[el]))
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
`PartialCombine(fids::Array{Int64,1},contprobs::Dict{Int64,Array{Float64,1}})`

Ok. Should take the output of contprobs and return the right output.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);;

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
`TotalCombine(D::DynState, location::Int64, comps::Dict{Int64,Int64}, contprobs::Dict{Int64,Array{Float64,1}})`
This will do what I want.  
Return all of the states, with the associated probabilities.  

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);;

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
nd = Dict{Int64,Int64}()
location = CompsDict(FindComps(dyn,dyn.all[1]), dyn, nd)
st_recs = Dict{Int64, NTuple{9,Int64}}()
StateRecord(nd, dyn, st_recs)

d1[1391330] = Dict{NTuple{10, Int64}, Float64}()
d1[1391330][TAddLevel(st_recs[1391330], 1)] = 0.0
d1[1391330][TAddLevel(st_recs[1391330], 2)] = 0.0
d1[1391330][TAddLevel(st_recs[1391330], 3)] = 0.0

cp1 = ContProbs(dyn.all[1].fid, st_recs, d1)

TotalCombine(dyn, 1, nd, cp1)


# try with two firms...


d1 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
d2 = Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  }()
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

d1[dyn.all[18].fid] = Dict{NTuple{10,Int64}, Float64}()
d1[dyn.all[18].fid][StateKey(dyn.all[18], 1)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 2)] = 0.0
d1[dyn.all[18].fid][StateKey(dyn.all[18], 3)] = 0.0
nd2 = Dict{Int64,Int64}()
location2 = CompsDict(FindComps(dyn, dyn.all[18]), dyn, nd2) # locations are 19 and 152
st_recs2 = Dict{Int64, NTuple{9,Int64}}()
StateRecord(nd2, dyn, st_recs2)

d1[dyn.all[19].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[19].fid][TAddLevel(st_recs2[dyn.all[19].fid], 1)] = 0.0
d1[dyn.all[19].fid][TAddLevel(st_recs2[dyn.all[19].fid], 2)] = 0.0
d1[dyn.all[19].fid][TAddLevel(st_recs2[dyn.all[19].fid], 3)] = 0.0

d1[dyn.all[152].fid] = Dict{NTuple{10, Int64}, Float64}()
d1[dyn.all[152].fid][TAddLevel(st_recs2[dyn.all[152].fid], 1)] = 0.0
d1[dyn.all[152].fid][TAddLevel(st_recs2[dyn.all[152].fid], 2)] = 1.0
d1[dyn.all[152].fid][TAddLevel(st_recs2[dyn.all[152].fid], 3)] = 2.0

cp2 = ContProbs(dyn.all[18].fid, st_recs2, d1)
dyn.all[19].nfids = [672285, 373510] # this correction should not be necessary.
compprobs = TotalCombine(dyn, 18, nd2, cp2)

"""

function TotalCombine(D::DynState,
                      location::Int64,
                      comps::Dict{Int64,Int64},# dict of fids/locations from nbs in ExactChoice.  
                      contprobs::Dict{Int64,Array{Float64,1}})
  # Create the outputs
  out05::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  push!(out05, ([0, 0, 0], 1.0))
  out515::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  push!(out515, ([0, 0, 0], 1.0))
  out1525::Array{Tuple{Array{Int64,1}, Float64}, 1} = Array{Tuple{Array{Int64,1}, Float64}, 1}()
  push!(out1525, ([0, 0, 0], 1.0))
  lastout::Dict{NTuple{9, Int64}, Float64} = Dict{NTuple{9, Int64}, Float64}()
  # Competitors' locations: should this call FindComps at all?  
  # comps::Array{Int64,1} = FindComps(D.all[location], D)
  # Measure the distances, apply GenStates to relevant segment.
  loc1 = [0.0 0.0 0.0]; 
  loc2 = [0.0 0.0 0.0];
  loc3 = [0.0 0.0 0.0];   
  for k in keys(comps)
    if comps[k] != location # this will ignore the firm and do only the neighbors.  
      k1 = comps[k] # access the location  
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
`StateRecord(neighbors::Dict{Int64,Int64},  D::DynState, outp::Dict{Int64,NTuple{9,Int64}})`
The idea would be that this would take a dict of fids and ints (locations) and then 
return a dict of fids/states, as StateRecord above does.
This would take the dict "locs" from ExactVal
This needs to catch the "special" value, I think.  That should be in "location".
  #### TESTS ####
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
Tex = EntireState(Array{Market,1}(), Dict{Int64, Market}(), Dict{Int64, Int64}());
CMakeIt(Tex, ProjectModule.fips);
FillState(Tex, ProjectModule.alldists, 50);
patients = NewPatients(Tex);

dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);;

test2 = Dict(3396057=>195, 3390720=>196 , 3396327=>197 , 3396189=>198, 2910645 => 11)
out_1 = Dict{Int64,NTuple{9,Int64}}()
StateRecord(test2, dyn,out_1)  
out_1 == Dict(3396057 => (0, 0, 1, 0, 0, 2, 1, 0, 0),3390720 => (0, 0, 0, 0, 1, 1, 1, 0, 1),3396327 => (0, 1, 0, 0, 0, 1, 1, 0, 1),3396189 => (0, 0, 0, 0, 1, 0, 1, 0, 2),2910645 => (0, 0, 0, 0, 0, 0, 0, 1, 3) ) # returns true 

v1 = FindComps(dyn, dyn.all[18])
push!(v1, 18)
d1 = Dict{Int64,Int64}()
test1 = Dict{Int64, NTuple{9,Int64}}()
CompsDict(v1, dyn, d1)
StateRecord(d1, dyn, test1)

test1 == Dict(670132 => (0, 0, 0, 1, 0, 0, 1, 0, 0),672285 => (0, 0, 0, 1, 0, 0, 0, 0, 0),373510 => (0, 0, 0, 0, 0, 0, 1, 0, 0)) # returns true. 


#Try with dyn.all[11]
v11 = FindComps(dyn, dyn.all[11])
push!(v11, 11)
d11 = Dict{Int64,Int64}()
test11 = Dict{Int64, NTuple{9,Int64}}()
CompsDict(v11, dyn, d11)
StateRecord(d11, dyn, test11)
test11 == Dict(3396057 => (0, 0, 1, 0, 0, 2, 1, 0, 0), 3390720 => (0, 0, 0, 0, 1, 1, 1, 0, 1), 3396327 => (0, 1, 0, 0, 0, 1, 1, 0, 1), 3396189 => (0, 0, 0, 0, 1, 0, 1, 0, 2),2910645 => (0, 0, 0, 0, 0, 0, 0, 1, 3))

"""
function StateRecord(neighbors::Dict{Int64,Int64},  D::DynState, outp::Dict{Int64,NTuple{9,Int64}})
  #outp::Dict{Int64, NTuple{9,Int64}} = Dict{Int64, NTuple{9,Int64}}() # records the neighbors of competing firms, relative to the firm we are computing EQ for, 
  for fid1 in keys(neighbors) # vector of locations of neighbors
    el1 = neighbors[fid1] # this is the location 
    intermed::Array{Int64,1} = zeros(Int64,9) 
    for fid2 in keys(neighbors) # check every other facility to get this state.  
      el2 = neighbors[fid2] # this has to allocate a tiny amount. 
      if el1!=el2
        d1::Float64 = distance(D.all[el1].lat, D.all[el1].long, D.all[el2].lat, D.all[el2].long)
          if (d1>0)&(d1<=5)
            if D.all[el2].level == 1
              intermed[1]+=1
            elseif D.all[el2].level == 2
              intermed[2]+=1
            elseif D.all[el2].level == 3
              intermed[3]+=1
            else
              #do nothing 
            end 
          elseif (d1>5)&(d1<=15)
            if D.all[el2].level == 1
              intermed[4]+=1
            elseif D.all[el2].level == 2
              intermed[5]+=1
            elseif D.all[el2].level == 3
              intermed[6]+=1
            else
              #do nothing 
            end 
          elseif (d1>15)&(d1<25)
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
  #return outp 
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


dyn = DynStateCreate(TexasEq, Tex, patients, ProjectModule.pcount);; 

p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

DSimNew(dyn.all[1].mk, dyn.all[1].fid, p1, p2)

PatientRev(dyn.all[1], p1, p2, 10) # recall that the scale factor is large so this should be small.
"""
function PatientRev(s::simh,
                    mpats::ProjectModule.patientcount,
                    ppats::ProjectModule.patientcount,
                    action::Int64)
    const scalefact::Float64 = 3.0e9,
    const alf1::Float64 = 8336.17,
    const alf2::Float64 = 36166.6,
    const alf3::Float64 = 16309.47,
    const gamma_1_385::Float64 = 20680.0, # ✓
    const gamma_2_385::Float64 = 42692.37, # ✓
    const gamma_3_385::Float64 = 20962.97, # ✓
    const gamma_1_386::Float64 = 81918.29, # X
    const gamma_2_386::Float64 = 74193.4, # X
    const gamma_3_386::Float64 = 99065.79, # X
    const gamma_1_387::Float64 = 30405.32, # X
    const gamma_2_387::Float64 = 49801.84, # X
    const gamma_3_387::Float64 = 22376.8, # X
    const gamma_1_388::Float64 = 10051.55, # ✓
    const gamma_2_388::Float64 = 19019.18, # X
    const gamma_3_388::Float64 = 33963.5, # X
    const gamma_1_389::Float64 = 29122.89, # X
    const gamma_2_389::Float64 = 14279.58, # X
    const gamma_3_389::Float64 = 20708.15, # X
    const gamma_1_390::Float64 = 22830.05, # X
    const gamma_2_390::Float64 = 6754.76, # X
    const gamma_3_390::Float64 = 3667.42, # ✓
    const gamma_1_391::Float64 = 9089.77, # X
    const gamma_2_391::Float64 = 8120.85, # X
    const gamma_3_391::Float64 = 1900.5, # ✓
    const mcaid385::Float64 = 151380.0,
    const mcaid386::Float64 = 48417.0,
    const mcaid387::Float64 = 18845.0,
    const mcaid388::Float64 = 7507.0,
    const mcaid389::Float64 = 9424.0,
    const mcaid390::Float64 = 4623.0,
    const mcaid391::Float64 = 3664.0 # to DRG mean added 3094 - avg reimbursement for DRGs 370-375 under TX Medicaid (2012)
    outp::Float64 = 0.0
    wtp::Float64 = FindWTP(s) # FIXME - can this give me a NaN? 
    if s.level == 1
      outp = alf1*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_1_385*(ppats.count385+mpats.count385) - gamma_1_386*(ppats.count386+mpats.count386) - gamma_1_387*(ppats.count387+mpats.count387) - gamma_1_388*(mpats.count388+ppats.count388) - gamma_1_389*(mpats.count389+ppats.count389) - gamma_1_390*(ppats.count390+mpats.count390) - gamma_1_391*(ppats.count391+mpats.count391)
    elseif s.level == 2
      outp = alf2*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_2_385*(ppats.count385+mpats.count385) - gamma_2_386*(ppats.count386+mpats.count386) - gamma_2_387*(ppats.count387+mpats.count387) - gamma_2_388*(mpats.count388+ppats.count388) - gamma_2_389*(mpats.count389+ppats.count389) - gamma_2_390*(ppats.count390+mpats.count390) - gamma_2_391*(ppats.count391+mpats.count391)
    elseif s.level == 3
      outp = alf3*wtp*(sum(ppats)) + mpats.count385*mcaid385 + mpats.count386*mcaid386 + mpats.count387*mcaid387 + mpats.count388*mcaid388 + mpats.count389*mcaid389 + mpats.count390*mcaid390 + mpats.count391*mcaid391 - gamma_3_385*(ppats.count385+mpats.count385) - gamma_3_386*(ppats.count386+mpats.count386) - gamma_3_387*(ppats.count387+mpats.count387) - gamma_3_388*(mpats.count388+ppats.count388) - gamma_3_389*(mpats.count389+ppats.count389) - gamma_3_390*(ppats.count390+mpats.count390) - gamma_3_391*(ppats.count391+mpats.count391)
    else # level = -999 (exit)
      outp = pi*scalefact  # FIXME - what is pi doing here?  
    end
    return outp/scalefact
end




"""
`StateKey(h::simh)`
Returns a tuple of the neighbors and level to use in `ExactVal()`

Note that this function should ONLY be used with the main firm, not with 
any neighbors, since this takes the full state directly from the DynState 
record.  For neighbors, use NStateKey. 
"""
function StateKey(h::simh, a::Int64)
  return (h.cns.level105, h.cns.level205, h.cns.level305, h.cns.level1515, h.cns.level2515, h.cns.level3515, h.cns.level11525, h.cns.level21525, h.cns.level31525, a)
end 

"""
`NStateKey(rstate::NTuple{9,Int64},lev::Int64)`

This function generates a state key for a firm which is a neighbor.  
Reasoning: Restricted states are generated as tuples.  StateKey takes 
full, unrestricted states from DynState, so only works for the main firm.
Then we use StateRecord to get restricted states as tuples.  Those should
form the state keys.  
For main firms, use StateKey.   
"""
function NStateKey(rstate::NTuple{9,Int64},lev::Int64)
  return (rstate..., lev)
end 



"""
`FindComps(D::DynState, arr::Array{Int64,1}, args::simh...)`
Finds the locations of the neighbors in the DynState.
Takes a variable number of simh.  
Appends to the input vector arr.
Returns a vector of the *locations* in D.all[] of the competitors of args...
Only pushes an address if it isn't already there.  
"""
function FindComps(D::DynState, arr::Array{Int64,1}, args::simh...)
  for (i, h) in enumerate(args)
    for el in h.nfids
      for loc in 1:size(D.all,1) 
        if D.all[loc].fid == el 
          if !in(loc, arr)
            push!(arr, loc)
          end 
        end  
      end 
    end 
  end 
end 

"""
`FindFids(D::Dynstate, locs::Array{Int64,1})`
Takes a list of simh and will return a list of all of the fids of their competitors.
This is just like FindComps but returns a different kind of object.  
"""
function FindFids(D::DynState, locs::Array{Int64,1})
  outp::Array{Int64,1} = Array{Int64,1}()
  for i in locs
    for el in D.all[i].nfids 
      push!(outp, el)
    end 
  end 
  return outp 
end 


"""
`MapCompState`

This is a fucking nightmare.  There must be a better way...
For EACH firm listed in ch, 
for EACH state in states,
for EACH fac in D.all[].mk.m.facs, 
IF fid in states 
update level according to states.
At the end update DUtil.  What the fuck...

TODO - finish this.

### Testing ### 

dyn = CounterObjects(50);
fids1 =  [1391330]
ch1 = [1]
states1 = [(1391330, 3)] # changing level

MapCompState(dyn, ch1, fids1, states1)

"""
function MapCompState(D::DynState, ch::Array{Int64,1}, fids::Array{Int64,1} , states::Array{Tuple{Int64,Int64}})
  for el in ch                     # these are locations in D.all 
    for zp in D.all[el].mk.m       # these are the zipcodes at each D.all[el]
      for f in zp.facs             # these are the facilities in the zip.
        if f.fid in(f, fids)       # update those which are relevant ONLY 
          GetNewLevel(f, states)   # calls this to fix the level and record that it must be updated.
        end 
      end 
    end 
  end 
end 

"""
`GetNewLevel(f::shortrec, states::Array{Tuple{Int64,Int64}})`
This is dumb.  Takes a row of (fid,level) tuples and a shortrec, and changes the level.
"""
function GetNewLevel(f::shortrec, states::Array{Tuple{Int64,Int64}})
  for el in states                       # Set of (fid,level) tuples
    if el[1] == f.fid 
      f.level = el[2]
      f.tbu = true                       # signals that it must be updated by UpdateDUtil.
    end 
  end 
end 



"""
`CompsDict(arr1::Array{Int64,1}, D::DynState, outp::Dict{Int64,Int64})`
Takes and array of ints which are locations in D.all and returns 
a dict{Fid, Location}.  For use with FindComps, which returns an array.
NB: does NOT append the location of the args... to FindComps.

### TESTING ###

testcd = Dict{Int64, Int64}()
CompsDict(FindComps(dyn, dyn.all[11]), dyn, testcd)
testcd == Dict(3396057 => 195, 3390720 => 196, 3396327 => 197, 3396189 => 198) # returns true.
"""
function CompsDict(arr1::Array{Int64,1}, D::DynState, outp::Dict{Int64,Int64})
  for el in arr1
    outp[D.all[el].fid] = el
  end
end 


"""
`UpdateD(h::simh)`
Updates the utility component.  The deterministic part.
This is used just for the main firm in the deterministic counterfactual.
This will *only* change by the amount of a switched level.  Does not recompute.
"""
function UpdateD(h::simh)
  if h.level == 1 # this item varies.
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
The change works like: UtilUp(c, fid, actual, current)
where "actual" is the field dyn.all[].actual - the immutable real level,
and "current" is the field dyn.all[].level - this can and may change.  
"""
function UtilUp(c::cpats, 
                fid::Int64, 
                actual::Int64,  # this one should be permanent.
                current::Int64) # this one should vary.
  const inteninter_med::Float64 = -0.572397,
  const interinten_med::Float64 = 0.572397,
  const inten_med::Float64 = 1.34994,
  const inter_med::Float64 = 0.777542,
  const inteninter_p::Float64 = 0.3197218,
  const interinten_p::Float64 = -0.3197218,
  const inten_p::Float64 = 1.18599,
  const inter_p::Float64 = 0.866268
  const indx_m::Int64 = findfirst(c.mutils[1,:], fid)
  const indx_p::Int64 = findfirst(c.putils[1,:], fid)
  if (actual == 1)&(current == 1)
    # do nothing.
  elseif (actual == 1)&(current == 2) # definitely want addition
    c.putils[2,indx_p] += inter_p
    c.mutils[2,indx_m] += inter_med
  elseif (actual == 1)&(current == 3) # definitely want addition
    c.putils[2,indx_p] += inten_p
    c.mutils[2,indx_m] += inten_med 
  elseif (actual == 2)&(current == 1) # This should be subtraction.
    c.putils[2,indx_p] -= inter_p
    c.mutils[2,indx_m] -= inter_med
  elseif (actual == 2)&(current == 2)
    # do nothing.
  elseif (actual == 2)&(current == 3) # definitely want adding here. 
    c.putils[2,indx_p] += inteninter_p 
    c.mutils[2,indx_m] += inteninter_med
  elseif (actual == 3)&(current == 1) # definitely want subtraction.
    c.putils[2,indx_p] -= inten_p 
    c.mutils[2,indx_m] -= inten_med 
  elseif (actual == 3)&(current == 2)
    c.putils[2,indx_p] += interinten_p   # FIXME - is this subtracted?
    c.mutils[2,indx_m] += interinten_med # FIXME - is this subtracted?
  elseif (actual == 3)&(current == 3)
    # do nothing.
  end 
end 


"""
`UtilDown(h::simh)`
Adjusts the deterministic utility back down.  
Now changed so that this function updates the level back to the actual value.  
"""
function UtilDown(h::simh)
  if h.level == 1
    if h.actual == 3
      for el in 1:size(h.mk.m,1) # this is an array of cpats
        UtilUp(h.mk.m[el], h.fid, 1, 3) # the syntax is: current level is 1, but the level is actually 3.
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
  # reset the level AFTER this is called:
  h.level = h.actual  
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
`DictCopy(d1::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, d2::Dict{ Int64, Dict{NTuple{10, Int64},  Float64}  }, alpha::Float64)`
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
DictClean(test1)
DictClean(test2)

test1[dyn.all[6].fid][(0,0,0,0,0,0,0,0,0,1)] = 20 #assign a value.

DictCopy(test2, test1, 0.5)

test2[dyn.all[6].fid][(0,0,0,0,0,0,0,0,0,1)] == 10.0 # should return true

"""
function DictCopy(d1::Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }, #permanent / outvals
                  d2::Dict{ Int64, Dict{NTuple{10, Int64},  Float64}  }, # temporary / tempvals
                  alpha::Float64)
  for k1 in keys(d1) # these are fids 
    for k2 in keys(d1[k1]) # these are neighbor state/level keys at the hospital level.  
      if !haskey(d2[k1],k2 ) # if the temporary doesn't have the key 
        # do nothing???  
      else # FIXME - there should be a test here concerning whether of these is zero, eg. the temp.  
        if d1[k1][k2] > 0 # don't copy states with zero value yet. - nope, this won't work.  Recall the initial value...  
          d1[k1][k2] = alpha*d2[k1][k2] + (1-alpha)*d1[k1][k2] # this should copy from the temp to the permanent.   
        end 
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
          outp[TupleSmash(i,j,k,1)] = 0.01 
        end
        if !haskey(outp, TupleSmash(i,j,k,2))
          outp[TupleSmash(i,j,k,2)] = 0.01
        end
        if !haskey(outp, TupleSmash(i,j,k,3))
          outp[TupleSmash(i,j,k,3)] = 0.01
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
`KeytoTuple(x::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64})`
The visit record is a giant tuple - this cuts that into (neighbors, level) and (action).
The former is a tuple, the latter a scalar.
"""
function KeytoTuple(x::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64})
  return x[1:end-1], x[end]
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

FIXME - here the shock should be added, I think.  
"""
function ChooseAction(h::simh)
  act::Int64 = sample(h.visited[KeyCreate(h.cns, h.level)].psi[1,:], Weights(h.visited[KeyCreate(h.cns, h.level)].psi[2,:]))
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









#=


=#
