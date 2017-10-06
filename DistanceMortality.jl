using DataFrames  
using CSV 


mutable struct DR # this is... how many patients traveled what distances from a zip. 
    p::patientcount 
    d::Float64 
end 


"""
`Finder`
Finds fid in d.all.
Returns 0 if not there.  
"""
function Finder(d::DynState, f::Int64)
  rr = 0
  for i = 1:size(d.all,1)
    if d.all[i].fid == f 
      rr = i 
      break  
    end 
  end
  return rr 
end 


"""
`DMMapCompState`

This is dumb.  Second version of this function just for this part to avoid breaking the first one.  
"""
function DMMapCompState(D::DynState, locs::Dict{Int64,Int64}, ch::Array{Int64,1}, fids::Array{Int64,1} , states::Array{Tuple{Int64,Int64}})
  if (length(states)>1)||(states[1][1] != 0)
    for el in ch                                                # these are locations in D.all - but there should be only one ALWAYS.
      for tp in states                                          # Levels need to be updated in the D - since these levels are drawn in UtilUp.
        if haskey(locs, tp[1])
          # do nothing 
        else 
          locs[tp[1]] = Finder(D, tp[1])
        end 
        if (D.all[locs[tp[1]]].level != tp[2])                  # level changes!
          for zp in D.all[el].mk.m                              # these are the zipcodes at each D.all[el]
            UtilUp(zp, tp[1], D.all[locs[tp[1]]].level, tp[2]; audit = false)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
          end  
          D.all[locs[tp[1]]].actual = D.all[locs[tp[1]]].level  # this is "storing where I was to fix it later"  
          D.all[locs[tp[1]]].level = tp[2]                      # this is "recording where I am" NB: timing of this line matters for utility update.  
        end 
      end 
      for zp in D.all[el].mk.m
        WTPNew(zp.putils, zp.pwtp)                              # once all utilities have been changed, fix the WTP.  
      end 
    end 
  end 
end 



"""
`MktDistance(d::DynState, 
                     chunk::Array{Int64,1}, 
                     all_locs::Dict{Int64,Int64}, 
                     conf::Array{Tuple{Int64,Int64},2}, 
                     medcounts::Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }, 
                     privcounts::Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}})`

Records all distances traveled to all firms in all zips around the firm in ch.
Records all patient volumes and choices for each firm in a vector of DR's, each 
containing the total number of patients choosing and total distances traveled
by that group.  Takes the argument `conf` which is a specific market configuration.
Also needs `all_locs` - fids and locations in `dyn.all[]` of the firms in `conf`.
Operates in place on medcounts and privcounts. 

## Testing ##

dyn = CounterObjects(1);

# Small Market 
medcounts = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }()
privcounts = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}()
chunk = [1];
conf = [(3490795, 1) (1391330, 2)]
MktDistance(dyn, [1], conf, medcounts, privcounts)

# Medium Market 
medcounts = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }()
privcounts = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}()
chunk = [11];
conf2 = [(3396057,1), (3390720,1), (3396327,1), (3396189,1), (2910645, 3)]
MktDistance(dyn, [11], conf2, medcounts, privcounts)

# Large Market 
medcounts1 = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }();
privcounts1 = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}();
chunk = [245];
conf2 = [(4530190,3) (4916068,1) (4916029,1) (4536048,1) (4530200,1) (4536337,1) (4530170,1) (4536338,1) (4536253,1)];
MktDistance(dyn, [245], conf2, medcounts1, privcounts1);

medcounts2 = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }();
privcounts2 = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}();
chunk = [245];
conf2 = [(4530190,1) (4916068,3) (4916029,3) (4536048,3) (4530200,3) (4536337,3) (4530170,3) (4536338,3) (4536253,3)];
MktDistance(dyn, [245], conf2, medcounts2, privcounts2);

# Make everyone exit:
medcounts3 = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }();
privcounts3 = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}();
chunk = [245];
conf2 = [(4530190,1) (4916068,999) (4916029,999) (4536048,999) (4530200,999) (4536337,999) (4530170,999) (4536338,999) (4536253,999)];
MktDistance(dyn, [245], conf2, medcounts3, privcounts3);

v1 = TakeAverage(dyn, medcounts1, privcounts1, 4530190)
v2 = TakeAverage(dyn, medcounts2, privcounts2, 4530190)
v3 = TakeAverage(dyn, medcounts3, privcounts3, 4530190)

  for i = 1:size(d.all, 1)
    if d.all[i].fid == conf[1][1]
        for j = 1:size(dyn.all[i].mk.m,1)
            for k = 1:size(dyn.all[i].mk.m[j].putils[1,:], 1)
                if dyn.all[i].mk.m[j].putils[1,k] == conf[1][1]
                    println(dyn.all[i].mk.m[j].zp, "  ", dyn.all[i].mk.m[j].putils[2,k])
                end 
            end 
        end 
    end 
  end
"""
function MktDistance(d::DynState, 
                     chunk::Array{Int64,1},  
                     conf::Array{Tuple{Int64,Int64}}, # this does take a configuration argument.  
                     medcounts::Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }, 
                     privcounts::Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}})
  k = d.all[chunk[1]].fid 
  all_locs::Dict{Int64,Int64} = Dict{Int64,Int64}()
  neighbors::Array{Int64,1} = Array{Int64,1}()
  FindComps(d, neighbors, d.all[chunk[1]])
  push!(neighbors, chunk[1])
  CompsDict(neighbors, d, all_locs)
  state = TotalFix(GiveState(d, chunk, all_locs, conf, d.all[all_locs[k]].cns), d, chunk)
  DMMapCompState(d, all_locs, chunk, FindFids(d, chunk), conf) # This can include a state for the firm in chunk.
  medcounts[state] =  Dict{Int64,Array{DR,1} }()
  privcounts[state] =  Dict{Int64,Array{DR,1} }()
  pcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  mcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  fds = Array{Int64,1}()
  temparr = zeros(2, 12)
  for el in conf 
    push!(fds, el[1])
  end 
  medcounts[state][d.all[all_locs[k]].fid] =  Array{DR,1}()
  privcounts[state][d.all[all_locs[k]].fid] =  Array{DR,1}()
  d1 = Dict{Int64,patientcount}()
  d2 = Dict{Int64,patientcount}()
  for i = 1:size(d.all[all_locs[k]].mk.m,1)
    # these compute the whole market demand for the state, i.e., the tuple.
    TotalMktDemand(d.all[chunk[1]].mk.m[i].putils, temparr, d1, PatExpByType(d.all[chunk[1]].mk.m[i].pcounts, true))
    TotalMktDemand(d.all[chunk[1]].mk.m[i].mutils, temparr, d2, PatExpByType(d.all[chunk[1]].mk.m[i].pcounts, false))
    for fr = 1:size(d.all[all_locs[k]].mk.m[i].putils[1,:],1) # copy the pcount and mcount for each firm. loop over firms !
      fdd = d.all[all_locs[k]].mk.m[i].putils[1,fr]
      if fdd != 0 # this skips the OO.
        mc1, pc1 = CopyCount(d1[fdd], d2[fdd]) 
        # now measure the distance.
        dis = DistanceGet(d, fdd, d.all[all_locs[k]].mk.m[i].lat, d.all[all_locs[k]].mk.m[i].long)
        # push the copied p and m 
        if !haskey(medcounts[state], fdd)
          medcounts[state][fdd] = Array{DR,1}()
        end  
        if !haskey(privcounts[state], fdd)
          privcounts[state][fdd] = Array{DR,1}()
        end 
        push!(medcounts[state][fdd], DR(mc1, dis))
        push!(privcounts[state][fdd], DR(pc1, dis))
        # reset patient counts 
        ResetP(pcount)
        ResetP(mcount) 
      end 
    end 
    CleanMktDemand(d1)
    CleanMktDemand(d2)
  end 
  ResetCompState(d, all_locs, chunk, FindFids(d, chunk), conf) # set it back 
end 



function MktProf(d::DynState, 
                 chunk::Array{Int64},  
                 conf::Array{Tuple{Int64,Int64}}, # this does take a configuration argument.  
                 medps::Dict{Int64, Array{Float64,1}} , 
                 privps::Dict{Int64, Array{Float64,1}})
  loc = Finder(d, chunk[1])
  k = d.all[loc].fid 
  all_locs::Dict{Int64,Int64} = Dict{Int64,Int64}()
  neighbors::Array{Int64,1} = Array{Int64,1}()
  FindComps(d, neighbors, d.all[loc])
  push!(neighbors, loc)
  CompsDict(neighbors, d, all_locs)
  state = TotalFix(GiveState(d, chunk, all_locs, conf, d.all[all_locs[k]].cns), d, chunk)
  DMMapCompState(d, all_locs, chunk, FindFids(d, chunk), conf) # This can include a state for the firm in chunk.
  
  
  fds = Array{Int64,1}()
  temparr = zeros(2, 12)
  for el in conf 
    push!(fds, el[1])
  end 


  d1 = Dict{Int64,patientcount}()
  d2 = Dict{Int64,patientcount}()
  for i = 1:size(d.all[all_locs[k]].mk.m,1)
    # these compute the whole market demand for the state, i.e., the tuple.
    TotalMktDemand(d.all[chunk[1]].mk.m[i].putils, temparr, d1, DrawAll(d.all[chunk[1]].mk.m[i].pcounts) )
    TotalMktDemand(d.all[chunk[1]].mk.m[i].mutils, temparr, d2, DrawAll(d.all[chunk[1]].mk.m[i].pcounts) )
    # Profit here.  

    CleanMktDemand(d1)
    CleanMktDemand(d2)
  end 
  ResetCompState(d, all_locs, chunk, FindFids(d, chunk), conf) # set it back 
end 





"""
`PopAvgDist`

Computes two average distances: one weighted by population, the other weighted by number of zip codes.

TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 1);
NewPatients(TexasEq);

48201 Houston, 48453 Travis, 48029 Bexar, 48141 El Paso, 48439 Tarrant, 48113 Dallas 

for cnty in [48201 48453 48029 48113 48439 48141]
  sm = 0.0
  for el in keys(TexasEq.mkts[cnty].collection) 
    sm += PopAvgDist(dyn, el)[1]
  end 
  println(cnty, "  ", sm/TexasEq.mkts[cnty].collection.count, "    ", TexasEq.mkts[cnty].collection.count)
end 
"""
function PopAvgDist(d::DynState, fid::Int64)
  k = Finder(d, fid)
  td = 0.0                   # total distance*total people 
  tp = 0.0                   # total people 
  tz = size(d.all[k].mk.m,1) # total zips 
  vtd = 0.0                  # unweighted distance  
  for i = 1:size(d.all[k].mk.m,1)
    sm = d.all[k].mk.m[i].pcounts.u385 + d.all[k].mk.m[i].pcounts.u386 + d.all[k].mk.m[i].pcounts.u387 + d.all[k].mk.m[i].pcounts.u388 + d.all[k].mk.m[i].pcounts.u389 + d.all[k].mk.m[i].pcounts.u390 + d.all[k].mk.m[i].pcounts.u391
    smm = d.all[k].mk.m[i].mcounts.u385 + d.all[k].mk.m[i].mcounts.u386 + d.all[k].mk.m[i].mcounts.u387 + d.all[k].mk.m[i].mcounts.u388 + d.all[k].mk.m[i].mcounts.u389 + d.all[k].mk.m[i].mcounts.u390 + d.all[k].mk.m[i].mcounts.u391
    ds1 = distance(d.all[k].mk.m[i].lat, d.all[k].mk.m[i].long, d.all[k].lat, d.all[k].long)
    td += ds1*(sm+smm)
    tp += (sm+smm)
    vtd += ds1 
  end 
  return round(td/tp,3), round(vtd/tz,3)
end





"""
`TotalFix`
Takes the state tuple and adds the own state to it.
This is a more intuitive way of representing the state of the whole market.  
"""
function TotalFix(t::NTuple{9,Int64}, dyn::DynState, ch::Array{Int64,1})
  if dyn.all[ch[1]].level == 1
    return (t[1]+1, t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9])
  elseif dyn.all[ch[1]].level == 2
    return (t[1], t[2]+1, t[3], t[4], t[5], t[6], t[7], t[8], t[9])
  elseif dyn.all[ch[1]].level == 3
    return (t[1], t[2], t[3]+1, t[4], t[5], t[6], t[7], t[8], t[9])
  end 
end 



"""
`TotalMktDemand`
Here I need to compute the demand for everyone in the zip actually...
Then I want to take this over all zips.
EACH ZIP must include...
For EACH hospital
A set of DR's. 
So it must include a dict of fid,patientcount.
d1 is {fid, Array{DR,1}}


dyn = CounterObjects(1);
chunk = [1];
p1 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
p2 = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)

all_locs = Dict(3490795=>1, 1391330=>90)


tempa = zeros(2,12)
d1 = Dict{Int64,patientcount}()
for el in dyn.all[10].mk.m[1].putils[1,:]
  d1[el] = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
end 
TotalMktDemand(dyn.all[120].mk.m[1].putils, tempa, d1, PatExpByType(dyn.all[120].mk.m[1].pcounts, true))

## Can test: ##
sum(PatExpByType(dyn.all[11].mk.m[1].pcounts, true))
tp = 0.0
for k1 in keys(d1)
  tp += sum(d1[k1])
end
This works in the cases I checked.  
"""
function TotalMktDemand(inparr::Array{Float64,2},temparr::Array{Float64,2}, d1::Dict{Int64, patientcount},pp::patientcount)
  WTPNew(inparr, temparr)                                           # updates temparr - now holds fids/WTP.
  for i = 1:size(temparr,2)                                         # this will loop over some zeros, but those will have zero mkt share.
    if !haskey(d1, temparr[1,i])
      d1[temparr[1,i]] = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)  # add the key if it happens to not be there.  
    end 
    d1[temparr[1,i]].count385 += pp.count385*temparr[2,i]           # recall that temparr will have a fraction in it.  
    d1[temparr[1,i]].count386 += pp.count386*temparr[2,i]
    d1[temparr[1,i]].count387 += pp.count387*temparr[2,i]
    d1[temparr[1,i]].count388 += pp.count388*temparr[2,i]
    d1[temparr[1,i]].count389 += pp.count389*temparr[2,i]
    d1[temparr[1,i]].count390 += pp.count390*temparr[2,i]
    d1[temparr[1,i]].count391 += pp.count391*temparr[2,i]
  end 
end 

"""
`CleanMktDemand`
Takes market demand and cleans it - sets all back to zero.  
"""
function CleanMktDemand(d1::Dict{Int64, patientcount})
  for k1 in keys(d1)
    d1[k1].count385 = 0.0
    d1[k1].count386 = 0.0
    d1[k1].count387 = 0.0
    d1[k1].count388 = 0.0
    d1[k1].count389 = 0.0
    d1[k1].count390 = 0.0
    d1[k1].count391 = 0.0
  end   
end 


"""
`DistanceGet`
Finds the location of f, the computes the distance, then returns a float with the distance.  
"""
function DistanceGet(d::DynState, f, lat::Float64, long::Float64)
  fid = 0
  dist = 0.0
  for i = 1:size(d.all, 1)
    if d.all[i].fid == f 
      fid = i 
      break 
    end 
  end 
  if fid != 0
    dist += distance(lat, long, d.all[fid].lat, d.all[fid].long)
  end 
  return dist 
end 


"""
`ResetP(pp::patientcount)`
Sets patientcount back to 0.
"""
function ResetP(pp::patientcount)
 pp.count385 = 0.0
 pp.count386 = 0.0
 pp.count387 = 0.0
 pp.count388 = 0.0
 pp.count389 = 0.0
 pp.count390 = 0.0
 pp.count391 = 0.0
end 

"""
`CopyCount(pc::patientcount, mc::patientcount)`
Copies the patientcounts.  Allocates new ones.  
"""
function CopyCount(pc::patientcount, mc::patientcount)
    p1::patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    p2::patientcount = patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    # First. 
    p1.count385 += pc.count385 
    p1.count386 += pc.count386 
    p1.count387 += pc.count387 
    p1.count388 += pc.count388 
    p1.count389 += pc.count389 
    p1.count390 += pc.count390 
    p1.count391 += pc.count391 
    # Second 
    p2.count385 += mc.count385 
    p2.count386 += mc.count386 
    p2.count387 += mc.count387 
    p2.count388 += mc.count388 
    p2.count389 += mc.count389 
    p2.count390 += mc.count390 
    p2.count391 += mc.count391 
    return p1, p2 
end 


"""
`TakeAverage(mc::Dict{NTuple{10,Int64}, Dict{Int64,Array{DR,1} } }, pc::Dict{NTuple{10,Int64}, Dict{Int64, Array{DR,1}}})`
Compute the average distances traveled per market config.  

dyn = CounterObjects(1);

medcounts = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }();
privcounts = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}();
conf2 = [(4530190,1) (4916068,3) (4916029,3) (4536048,3) (4530200,3) (4536337,3) (4530170,3) (4536338,3) (4536253,3)];
MktDistance(dyn, [245], conf2, medcounts, privcounts);

TakeAverage(dyn, medcounts, privcounts, 4530190)

Change this: return BOTH averages to that one facility, but also over the whole market.  

"""
function TakeAverage(d::DynState, mc::Dict, pc::Dict, f::Int64)
  loc = 0
  for i = 1:size(d.all,1)
    if d.all[i].fid == f 
      loc = i 
    end 
  end   
  cols = 13                                    # outp... fid, nine states, one level, one count of patients, one average distance
  rows = size(d.all[loc].nfids,1)+2            # neighbors, main firm, total for all.    
  outp = zeros(rows, cols)
  rc = 1                                       # row counter
  for k1 in keys(mc)                           # the state - i.e., the market configuration 
    for k2 in keys(mc[k1])                     # the fids.
      pats::Float64 = 0.0
      ds::Float64 = 0.0
      if (in(k2, d.all[loc].nfids))||(k2 == f) # checks to make sure the firm is in the market.
        for i = 1:size(mc[k1][k2],1)
          a1, b1 = DREX(mc[k1][k2][i])
          pats += a1                           # records total number of medicare patients.
          ds += a1*b1                          # records distance traveled by patients in that zip.
        end 
        for i = 1:size(pc[k1][k2],1)
          a1, b1 = DREX(pc[k1][k2][i])
          pats += a1                           # records total number of privately insured patients.
          ds += a1*b1                          # records distance traveled by patients in that zip.
        end
        outp[rc,1] = k2                        # writes out which FID is being checked.  
        for j = 1:length(k1)
          outp[rc, j+1] = k1[j]                # this records the state (a 9 tuple) as nine columns  
        end 
        # outp[rc,11] =                        # Level can be omitted, but then one column is always zero.
        outp[rc,12] += pats 
        outp[rc,13] += (ds/max(pats,1))        # this is: total number of miles traveled divided by all patients to that hospital.  
        rc += 1  
      end 
    end 
  end 
  totp = sum(outp[:,12])                       # now values over the whole market: how many patients over how many miles....
  sm_miles = 0.0
  for i = 1:size(outp,1)
    sm_miles += outp[i,12]*outp[i,13]
  end 
  outp[rows, 1] = 99999999                     # special fid 
  for j = 1:9
    outp[rows, j+1] = outp[rows-1, j+1]
  end 
  outp[rows,12] = totp 
  outp[rows,13] = sm_miles/totp 
  out1 = convert(DataFrame, outp)
  names!(out1, [:fid, :lev1_05, :lev2_05, :lev3_05, :lev1_515, :lev2_515, :lev3_515, :lev1_1525, :lev2_1525, :lev3_1525, :level, :patients, :avg_distance])
  return out1
end



"""
`Mortality(mc::Dict, pc::Dict, conf::Array{Tuple{Int64,Int64},2};
           mp_lin::Float64 = 0.005,
           mp_quad::Float64 = 0.0,
           nicad::Float64 = 0.06,
           fvlbw::Float64 = 0.014)`
Computes Medicaid and Private patient mortality among all hospitals.


## Testing ## 
dyn = CounterObjects(1);

medcounts2 = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }()
privcounts2 = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}()
chunk = [245];
conf2 = [(4530190,3), (4916068,3), (4916029,3), (4536048,3), (4530200,3), (4536337,3), (4530170,3), (4536338,3), (4536253,3)]
MktDistance(dyn, [245], conf2, medcounts2, privcounts2)

Mortality(medcounts2, privcounts2, conf2)

"""
function Mortality(mc::Dict, pc::Dict, conf::Array{Tuple{Int64,Int64}};
                   nicad::Float64 = 0.06,
                   fvlbw::Float64 = 0.014,
                   regionalize::Bool = false,
                   sp_fid::Int64 = 99999999)
  Ns::Int64 = 100                            # draws of mortality rate.  
  fidloc = 1                                 # fid  location 
  birthloc = 2                               # total births location 
  niculoc = 3                                # nicu admits location 
  vlbwloc = 4                                # total vlbw location 
  mploc = 5                                  # mortality prob location 
  mortloc = 6                                # total mortality  location 
  msdloc = 7                                 # mortality sd location 
  hhiloc = 8                                 # HHI  
  fds = Array{Int64,1}()
  rws = length(conf)+1
  for i = 1:(rws-1) 
    push!(fds, conf[i][1])                   # collect the fids. 
  end 
  cls = 8                                    # fid, number of admits, nicu admits, vlbw, mortality rate, mean total deaths, sd total deaths, hhi.
  outp::Array{Float64,2} = zeros(rws, cls)
  rc = 1                                     # counts rows!  in outp.  
  for k1 in keys(mc)                         # this is the market state 
    for k2 in keys(mc[k1])                   # these are the firms.
      if in(k2, fds)                         # these are the relevant firms.
        dths = zeros(Ns)                     # will hold the Ns draws of the mortality rate.  
        ct = 0.0                             # count of patients                         
        for j = 1:size(mc[k1][k2], 1)        # medicaid patients 
          a, b = DREX(mc[k1][k2][j])
          ct += a 
        end 
        for j = 1:size(pc[k1][k2], 1)        # private patients 
          a, b = DREX(pc[k1][k2][j])
          ct += a 
        end 
        cna = ct*nicad                       # count of nicu admitted patients.
        cvln = fvlbw*ct                      # count of vlbw patients 
        mps = 0.0 
        for m = 1:Ns
          mp = MortProb(cvln)
          mps += mp                          # track total of probs.
          if cvln*mp < 0.0
            println("prob: ", mp, " count ", cvln)
          end 
          dths[m] = cvln*mp
        end 
        outp[rc, fidloc] = k2                 # firm fid 
        outp[rc, birthloc] = ct               # birth count 
        outp[rc, niculoc] = cna               # nicu admits 
        outp[rc, vlbwloc] = cvln              # vlbw 
        outp[rc, mploc] = mps/Ns              # mort prob mean, over Ns draws.  
        outp[rc, mortloc] = mean(dths)        # total mort mean 
        outp[rc, msdloc] = std(dths)          # total mortality st. d.
        rc += 1
      end 
    end 
  end 
  if regionalize 
    dths = zeros(Ns)
    outp[rws, fidloc] = sp_fid                            # special fid in regionalized case 
    bc = sum(outp[:,birthloc])                            # sum total births
    outp[rws, birthloc] =                                 # record birth count 
    nic_ad = sum(outp[:,niculoc])                         # total nicu admits in market
    outp[rws, niculoc] = nic_ad          
    cvln = sum(outp[:,vlbwloc])                           # total vlbw in market
    outp[rws, vlbwloc] = cvln  
    mps = 0.0
    for m = 1:Ns
      mp = MortProb(cvln)
      mps += mp 
      dths[m] = cvln*mp
    end          
    outp[rws, mploc] = mps/Ns   
    outp[rws, mortloc] = mean(dths)                       # total mortality in market 
    outp[rws, msdloc] = std(dths)                         # standard deviation of market deaths. 
    hhi = 0.0
    for i = 1:(rws-1)
      hhi += (outp[i,birthloc]/bc)^2
    end 
    outp[rws, hhiloc] = hhi 
  else 
    dths = zeros(Ns)
    outp[rws, fidloc] = sp_fid                            # here fid will be 999999 - not regionalizing. 
    bc = sum(outp[:,2])                                   # compute total birth count
    outp[rws, birthloc] = bc                              # write birth count.
    outp[rws, niculoc] = sum(outp[:,3])                   # total nicu admits 
    outp[rws, vlbwloc] = sum(outp[:,4])                   # total vlbw 
    outp[rws, mploc] = mean(outp[1:(rws-1),5])            # mean mortality rate over all facilities.
    outp[rws, mortloc] = sum(outp[:,6])                   # total mean mortality 
    outp[rws, msdloc] = 0.0                               # not computing mean or sd of mortality rates.
    hhi = 0.0
    for i = 1:(rws-1)
      hhi+=(outp[i,birthloc]/bc)^2
    end 
    outp[rws, hhiloc] = hhi  
  end 
  out1 = convert(DataFrame, outp)                    # returning a dataframe just to use the column naming capability
  names!(out1, [:fid, :totalbirths, :nicu_admits, :vlbw, :mean_mort_rate, :mean_mortality, :std_mortality, :hhi])
  return out1
end 


"""
`MortProb(v;lp::Float64 = 0.5, qp::Float64 = 0.001)`
Returns a volume-implied hospital specific mortality rate.
Uses estimates from the TX Birth Certificate Data.

"""
function MortProb(v::T) where T <: Real
  lin_μ = 0.0
  lin_σ = 0.0
  if (v >= 0.0) & (v<20.0)           # these are marginal effects at mean values from data.  
    lin_μ = 0.151
    lin_σ = 0.0005
  elseif (v>=20.0) & (v <40.0)
    lin_μ = 0.0145
    lin_σ = 0.0004
  elseif (v>=40.0) & (v <60.0)
    lin_μ = 0.0138
    lin_σ = 0.0003
  elseif (v>=60.0) & (v <80.0)
    lin_μ = 0.0132
    lin_σ = 0.0003
  elseif (v>=80.0) & (v <100.0)
    lin_μ = 0.0126
    lin_σ = 0.0002
  elseif (v>=100.0) & (v <120.0)
    lin_μ = 0.0121
    lin_σ = 0.0002
  elseif (v>=120.0) & (v <140.0)
    lin_μ = 0.0115
    lin_σ = 0.0002
  elseif (v>=140.0) & (v <160.0)
    lin_μ = 0.0110
    lin_σ = 0.0002
  elseif (v>=160.0) & (v <180.0)
    lin_μ = 0.0105
    lin_σ = 0.0002
  elseif (v>=180.0) & (v <200.0)
    lin_μ = 0.0100
    lin_σ = 0.0003
  elseif (v>=200.0) & (v <220.0)
    lin_μ = 0.0096
    lin_σ = 0.0003
  elseif (v>=220.0) & (v <240.0)
    lin_μ = 0.0091
    lin_σ = 0.0003
  elseif (v>=240.0) & (v <260.0)
    lin_μ = 0.0087
    lin_σ = 0.0003
  elseif (v>=260.0) & (v <280.0)
    lin_μ = 0.0083
    lin_σ = 0.0004
  elseif (v>=280.0) & (v <300.0)
    lin_μ = 0.0079
    lin_σ = 0.0004
  elseif (v>=300.0) & (v <320.0)
    lin_μ = 0.0075
    lin_σ = 0.0004
  elseif (v>=320.0) & (v <340.0)
    lin_μ = 0.0072
    lin_σ = 0.0005
  elseif (v>=340.0) & (v <360.0)
    lin_μ = 0.0068
    lin_σ = 0.0005
  elseif (v>=360.0) & (v <380.0)
    lin_μ = 0.0065
    lin_σ = 0.0005
  elseif (v>=380.0) & (v <400.0)
    lin_μ = 0.0062
    lin_σ = 0.0005
  elseif (v>=400.0) & (v <420.0)
    lin_μ = 0.0059
    lin_σ = 0.0005
  elseif (v>=420.0) & (v <440.0)
    lin_μ = 0.0056
    lin_σ = 0.0005
  elseif (v>=440.0) & (v <460.0)
    lin_μ = 0.0053
    lin_σ = 0.0005
  else # v > 460   
    lin_μ = 0.0051
    lin_σ = 0.0005
  end 
  lind = Distributions.Normal(lin_μ, lin_σ)
  lp = Distributions.rand(lind)             # draw from the distribution 
  if lp < 0.0
    lp = abs(lp)                            # This is messing with the distribution a bit, but this should rarely be negative.
  end 
  return lp                                 # quadratic component could be added.
end 





"""
`DREX(d::DR)`
Returns sum of patients and distance.
"""
function DREX(d::DR)
  return sum(d.p), d.d 
end 


"""
`FindThem`
dyn = CounterObjects(1);
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 1);
NewPatients(TexasEq);
FindThem(dyn, TexasEq)

Do this by county.  55 counties have no neighbors.
99 counties have neighbors.
- If neighbor, get fids, call MktDistance
- Design config for that.  care about mortality only.  
- Do as given, then with one lev 3 
- compare 

Big Markets:
Austin 4536253

"""
function FindThem(d::DynState, es::EntireState)
  # Need an output holder for the whole set of outputs.  
  outdf = DataFrame( fid = @data([0.0]), 
                     lev1_05= @data([0.0]), 
                     lev2_05= @data([0.0]), 
                     lev3_05= @data([0.0]), 
                     lev1_515= @data([0.0]), 
                     lev2_515= @data([0.0]), 
                     lev3_515= @data([0.0]), 
                     lev1_1525= @data([0.0]), 
                     lev2_1525= @data([0.0]), 
                     lev3_1525= @data([0.0]), 
                     level = @data([0.0]), 
                     patients = @data([0.0]), 
                     avg_distance = @data([0.0]), 
                     totalbirths = @data([0.0]), 
                     nicu_admits = @data([0.0]), 
                     vlbw = @data([0.0]), 
                     mean_mort_rate = @data([0.0]), 
                     mean_mortality = @data([0.0]), 
                     std_mortality = @data([0.0]),
                     hhi=@data([0.0]), 
                     fipscode = @data([0.0]), 
                     counterfactual = @data([0.0])) 
  # Do these once - will do markets with more than neighbor at a firm.  
  todo = Dict{Int64, Bool}()
  for i in keys(es.mkts)
    for j in keys(es.mkts[i].collection)
      todo[i] = false 
      if (length(es.mkts[i].collection[j].hood) >0)&!(todo[i])             # this will do only markets where firms have multiple neighbors
        todo[i] = true 
      end 
    end 
  end 
  for k1 in keys(todo)                                                     # fipscodes for places to do 
    actual_arr::Array{Tuple{Int64,Int64},1}=Array{Tuple{Int64,Int64},1}()  # actual configuration 
    new_arr::Array{Tuple{Int64,Int64},1}=Array{Tuple{Int64,Int64},1}()     # assigning one level 3
    fc = 0                                                                 # count facilities 
    medcts = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1}}}()              # dict for distance count 
    privcts = Dict{NTuple{9,Int64},Dict{Int64,Array{DR,1}}}()              # dict for distance count 
    for k2 in keys(es.mkts[k1].collection)
      ix = 0
      if todo[k1] 
        for k3 in keys(es.mkts[k1].collection)                             # keys twice...?
          if fc == 0
            push!(actual_arr, (es.mkts[k1].collection[k3].fid, es.mkts[k1].collection[k3].level))
            push!(new_arr, (es.mkts[k1].collection[k3].fid, 3))
            ix = Finder(d, k3)
            fc+=1
          else 
            push!(actual_arr, (es.mkts[k1].collection[k3].fid, es.mkts[k1].collection[k3].level))
            push!(new_arr, (es.mkts[k1].collection[k3].fid, 1))
          end 
        end 
        # Equilibrium Arrangement. 
        medcts = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1}}}()              # dict for distance count 
        privcts = Dict{NTuple{9,Int64},Dict{Int64,Array{DR,1}}}()              # dict for distance count
        MktDistance(d, [ix], actual_arr, medcts, privcts)                      # distances computed.
        df11 = join(TakeAverage(d, medcts, privcts, actual_arr[1][1]), Mortality(medcts, privcts, actual_arr), on = :fid, kind = :outer)
        county = zeros(size(df11,1), 3)                                        # this only needs to be done once.  
        for (idx,v) in enumerate(df11[:fid])
          county[idx, 3] = 1                                                   # signals what counterfactual.  
          county[idx, 2] = k1
          county[idx, 1] = v
        end 
        county = convert(DataFrame, county)
        names!(county, [:fid, :fipscode, :counterfactual])
        df21 = join(df11, county, on = :fid, kind = :outer)
        append!(outdf, df21)
        CleanDistDict(medcts)
        CleanDistDict(privcts)
        # Single Level 3
        medcts = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1}}}()              # dict for distance count 
        privcts = Dict{NTuple{9,Int64},Dict{Int64,Array{DR,1}}}()              # dict for distance count
        MktDistance(d, [ix], new_arr, medcts, privcts)                         # distances computed.
        df21 = join(TakeAverage(d, medcts, privcts, new_arr[1][1]), Mortality(medcts, privcts, new_arr), on = :fid, kind = :outer)
        for (idx,v) in enumerate(df11[:fid])
          county[idx, 3] = 2                                                   # signals what counterfactual.  
        end 
        df22 = join(df21, county, on = :fid, kind = :outer)
        append!(outdf, df22)
        CleanDistDict(medcts)
        CleanDistDict(privcts)
        # Single Regionalized Level 3
        medcts = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1}}}()              # dict for distance count 
        privcts = Dict{NTuple{9,Int64},Dict{Int64,Array{DR,1}}}()              # dict for distance count
        MktDistance(d, [ix], new_arr, medcts, privcts)                         # distances computed.
        df31 = join(TakeAverage(d, medcts, privcts, new_arr[1][1]), Mortality(medcts, privcts, new_arr; regionalize = true), on = :fid, kind = :outer)
        for (idx,v) in enumerate(df11[:fid])
          county[idx, 3] = 3                                                   # signals what counterfactual.  
        end 
        df23 = join(df31, county, on = :fid, kind = :outer)
        append!(outdf, df23)
        CleanDistDict(medcts)
        CleanDistDict(privcts)
        todo[k1] = false                                                       # once finished, don't do again.  
      end 
    end
  end 
  return  outdf
end 




"""


mutable struct DR # this is... how many patients traveled what distances from a zip. 
    p::patientcount 
    d::Float64 
end
"""
function CleanDistDict(d1::Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1}}})
  for k1 in keys(d1)
    for k2 in keys(d1[k1])
      for i = 1:size(d1[k1][k2],1)
        d1[k1][k2][i].p.count385 = 0.0
        d1[k1][k2][i].p.count386 = 0.0
        d1[k1][k2][i].p.count387 = 0.0
        d1[k1][k2][i].p.count388 = 0.0
        d1[k1][k2][i].p.count389 = 0.0
        d1[k1][k2][i].p.count390 = 0.0
        d1[k1][k2][i].p.count391 = 0.0
        d1[k1][k2][i].d = 0.0
      end 
    end 
  end 
end

#=
# Running this and saving:
using ProjectModule
include("DistanceMortality.jl")
dyn = CounterObjects(1);
TexasEq = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 1);
NewPatients(TexasEq);
df1 = FindThem(dyn, TexasEq)

# dealing with missing vals is not working...
for i = 1:size(df1,1)
  for j = 1:size(df1,2)
    if isna(df1[i,j])
      df1[i,j] = -1.0
    end
  end
end


CSV.write("reg_counter_results.csv", ab; header = true)

=#
