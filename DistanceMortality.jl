using DataFrames  



mutable struct DR # this is... how many patients traveled what distances from a zip. 
    p::patientcount 
    d::Float64 
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
conf2 = [(3396057,1) (3390720,1) (3396327,1) (3396189,1) (2910645, 3)]
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

v1 = TakeAverage(dyn, medcounts1, privcounts1, 4530190)
v2 = TakeAverage(dyn, medcounts2, privcounts2, 4530190)

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
                     conf::Array{Tuple{Int64,Int64},2}, # this does take a configuration argument.  
                     medcounts::Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }, 
                     privcounts::Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}})
  k = d.all[chunk[1]].fid 
  all_locs::Dict{Int64,Int64} = Dict{Int64,Int64}()
  neighbors::Array{Int64,1} = Array{Int64,1}()
  FindComps(d, neighbors, d.all[chunk[1]])
  push!(neighbors, chunk[1])
  CompsDict(neighbors, d, all_locs)
  state = TotalFix(GiveState(d, chunk, all_locs, conf, d.all[all_locs[k]].cns), d, chunk)
  MapCompState(d, all_locs, chunk, FindFids(d, chunk), conf) # This can include a state for the firm in chunk.
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
    d1[temparr[1,i]].count385 += pp.count385*temparr[2,i]
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
    dist += distance(lat, long, dyn.all[fid].lat, dyn.all[fid].long)
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


TakeAverage(dyn, medcounts, privcounts, 3490795)

"""
function TakeAverage(d::DynState, mc::Dict, pc::Dict, f::Int64)
  # this can be changed... only check those facilities which are actual neighbors of f.  
  # Otherwise distances include many other firms.  
  loc = 0
  for i = 1:size(d.all,1)
    if d.all[i].fid == f 
      loc = i 
    end 
  end   
  cols = 13                             # outp... fid, nine states, one level, one count of patients, one average distance
  rows = mc.count 
  outp = zeros(rows, cols)
  rc = 1                                # row counter
  for k1 in keys(mc)                    # the state 
    pats::Float64 = 0.0
    ds::Float64 = 0.0
    for k2 in keys(mc[k1])              # the other firms.
      if (in(d.all[loc].nfids, k2))||(k2 == f) # checks to make sure the firm is in the market.
        println(k2)
        for i = 1:size(mc[k1][k2],1)
          a1, b1 = DREX(mc[k1][k2][i])
          pats += a1                      # records total number of medicare patients.
          ds += a1*b1                     # records distance traveled by patients in that zip.
        end 
        for i = 1:size(pc[k1][k2],1)
          a1, b1 = DREX(pc[k1][k2][i])
          pats += a1                      # records total number of privately insured patients.
          ds += a1*b1                     # records distance traveled by patients in that zip.
        end
      end 
    end 
    # now put it in the output in some way... 
    outp[rc,1] = f 
    for j = 1:length(k1)
      outp[rc, j+1] = k1[j]
    end 
    # outp[rc,11] =  # Level can be omitted, but then one column is always zero.
    outp[rc,12] += pats 
    outp[rc,13] += (ds/max(pats,1))       # this is: total number of miles traveled divided by all patients to that hospital.  
    rc += 1
  end 
  return outp
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
conf2 = [(4530190,3) (4916068,3) (4916029,3) (4536048,3) (4530200,3) (4536337,3) (4530170,3) (4536338,3) (4536253,3)]
MktDistance(dyn, [245], conf2, medcounts2, privcounts2)

Mortality(medcounts2, privcounts2, conf2)

"""
function Mortality(mc::Dict, pc::Dict, conf::Array{Tuple{Int64,Int64},2};
                   nicad::Float64 = 0.06,
                   fvlbw::Float64 = 0.014,
                   regionalize::Bool = false,
                   sp_fid::Int64 = 99999999)
  Ns::Int64 = 100                            # draws of mortality rate.  
  fds = Array{Int64,1}()
  rws = length(conf)+1
  for i = 1:(rws-1) 
    push!(fds, conf[i][1])                   # collect the fids. 
  end 
  cls = 7                                    # fid, number of admits, nicu admits, vlbw, mortality rate, mean total deaths, sd total deaths.
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
        outp[rc, 1] = k2                     # firm fid 
        outp[rc, 2] = ct                     # birth count 
        outp[rc, 3] = cna                    # nicu admits 
        outp[rc, 4] = cvln                   # vlbw 
        outp[rc, 5] = mps/Ns                 # mort prob mean, over Ns draws.  
        outp[rc, 6] = mean(dths)             # total mort mean 
        outp[rc, 7] = std(dths)              # total mortality st. d.
        rc += 1
      end 
    end 
  end 
  if regionalize 
    dths = zeros(Ns)
    outp[rws, 1] = sp_fid                            # special fid in regionalized case 
    outp[rws, 2] = sum(outp[:,2])                    # total birth count 
    nic_ad = sum(outp[:,3])                          # total nicu admits in market
    outp[rws, 3] = nic_ad          
    vlbw_t = sum(outp[:,4])                          # total vlbw in market
    outp[rws, 4] = vlbw_t  
    mps = 0.0
    for m = 1:Ns
      mp = MortProb(cvln)
      mps += mp 
      dths[m] = cvln*mp
    end          
    outp[rws, 5] = mps/Ns   
    outp[rws, 6] = mean(dths)                        # total mortality in market 
    outp[rws, 7] = std(dths)                         # standard deviation of market deaths.  
  else 
    dths = zeros(Ns)
    outp[rws, 1] = sp_fid                            # here fid will be 999999 - not regionalizing.  
    outp[rws, 2] = sum(outp[:,2])                    # total birth count 
    outp[rws, 3] = sum(outp[:,3])                    # total nicu admits 
    outp[rws, 4] = sum(outp[:,4])                    # total vlbw 
    outp[rws, 5] = mean(outp[1:(rws-1),5])           # mean mortality rate over all facilities.
    outp[rws, 6] = sum(outp[:,6])/(rws-1)            # total mean mortality 
    outp[rws, 7] = 0.0                               # not computing mean or sd of mortality rates.  
  end 
  out1 = convert(DataFrame, outp)                    # returning a dataframe just to use the column naming capability
  names!(out1, [:fid, :totalbirths, :nicu_admits, :vlbw, :mean_mort_rate, :mean_mortality, :std_mortality])
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
  if v >= 0.0 & v<20.0           # these are marginal effects at mean values from data.  
    lin_μ = 0.151
    lin_σ = 0.0005
  elseif v>=20.0 & v <40.0
    lin_μ = 0.0145
    lin_σ = 0.0004
  elseif v>=40.0 & v <60.0
    lin_μ = 0.0138
    lin_σ = 0.0003
  elseif v>=60.0 & v <80.0
    lin_μ = 0.0132
    lin_σ = 0.0003
  elseif v>=80.0 & v <100.0
    lin_μ = 0.0126
    lin_σ = 0.0002
  elseif v>=100.0 & v <120.0
    lin_μ = 0.0121
    lin_σ = 0.0002
  elseif v>=120.0 & v <140.0
    lin_μ = 0.0115
    lin_σ = 0.0002
  elseif v>=140.0 & v <160.0
    lin_μ = 0.0110
    lin_σ = 0.0002
  elseif v>=160.0 & v <180.0
    lin_μ = 0.0105
    lin_σ = 0.0002
  elseif v>=180.0 & v <200.0
    lin_μ = 0.0100
    lin_σ = 0.0003
  elseif v>=200.0 & v <220.0
    lin_μ = 0.0096
    lin_σ = 0.0003
  elseif v>=220.0 & v <240.0
    lin_μ = 0.0091
    lin_σ = 0.0003
  elseif v>=240.0 & v <260.0
    lin_μ = 0.0087
    lin_σ = 0.0003
  elseif v>=260.0 & v <280.0
    lin_μ = 0.0083
    lin_σ = 0.0004
  elseif v>=280.0 & v <300.0
    lin_μ = 0.0079
    lin_σ = 0.0004
  elseif v>=300.0 & v <320.0
    lin_μ = 0.0075
    lin_σ = 0.0004
  elseif v>=320.0 & v <340.0
    lin_μ = 0.0072
    lin_σ = 0.0005
  elseif v>=340.0 & v <360.0
    lin_μ = 0.0068
    lin_σ = 0.0005
  elseif v>=360.0 & v <380.0
    lin_μ = 0.0065
    lin_σ = 0.0005
  elseif v>=380.0 & v <400.0
    lin_μ = 0.0062
    lin_σ = 0.0005
  elseif v>=400.0 & v <420.0
    lin_μ = 0.0059
    lin_σ = 0.0005
  elseif v>=420.0 & v <440.0
    lin_μ = 0.0056
    lin_σ = 0.0005
  elseif v>=440.0 & v <460.0
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
  return (lp*v)                             # quadratic component could be added.
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

Do this by county.  55 counties have no neighbors.
99 counties have neighbors.
- If neighbor, get fids, call MktDistance
- Design config for that.  care about mortality only.  
- Do as given, then with one lev 3 
- compare 

for i in keys(TexasEq.mkts)
  for j = 1:size(TexasEq.mkts[i].config,1)
    if length(TexasEq.mkts[i].config[j].hood)>7
      println(i, "  ", TexasEq.mkts[i].config[j].fid,"  ", length(TexasEq.mkts[i].config[j].hood))
    end
  end
end
ct = 0
for i in keys(TexasEq.mkts)
  for j = 1:size(TexasEq.mkts[i].config,1)
    if length(TexasEq.mkts[i].config[j].hood)==0
      println(i, "  ", TexasEq.mkts[i].config[j].fid,"  ", length(TexasEq.mkts[i].config[j].hood))
      ct+=1
    end
  end
end


Big Markets:
Austin 4536253

"""
function FindThem(d::DynState, es::EntireState)
  # Need an output holder for the whole set of outputs. 
  # TODO - get these dims right.  
  # cols: fipscode, fid, ... whatever's in results... , 
  rows = 1000
  cols = 20
  outp::Array{Float64,2}(rows, cols)
  # Do these once - will do markets with more than neighbor at a firm.  
  todo = Dict{Int64, Bool}()
  for i in keys(es.mkts)
    for j in keys(es.mkts[i].collection)
      todo[i] = false 
      if (length(es.mkts[i].collection[j].hood) >0)&!(todo[i])
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
    ix = 0
    for k2 in keys(es.mkts[k1].collection)
      if todo[k2]
        if fc == 0
          push!(actual_arr, (es.mkts[k1].collection[k2].fid, es.mkts[k1].collection[k2].level))
          push!(new_arr, (es.mkts[k1].collection[k2].fid, 3))
          for i = 1:size(d.all,1)
            if d.all[i].fid == k2 
              ix = i 
            end 
          end 
          fc+=1
        else 
          push!(actual_arr, (es.mkts[k1].collection[k2].fid, es.mkts[k1].collection[k2].level))
          push!(new_arr, (es.mkts[k1].collection[k2].fid, 1))
        end 
      end 
    end
    MktDistance(d, [ix], actual_arr, medcts, privcts) # distances computed.  
    Mortality(medcts, privcts, actual_arr)
    CleanDicts(medcts, privcts)

  end  
end 




#=


=#
