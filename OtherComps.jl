#Other paper computations: 

#=

Miscellaneous model computations for the paper:

=#



### Some work with types and different processes...

# l1 = addprocs()

# @everywhere mutable struct smth 
#     f::Array{Int64,1}
# end 

# @everywhere function adder(n::smth, i::Int64)
#     for j = 1:i 
#         push!(n.f, j)
#     end 
#     return n # this function needs an explicit return statement, otherwise remotecall_fetch won't return anything.
# end 

# @everywhere function getter(n::smth)
#     return n.f 
# end 


# @everywhere function buildsm()
# return smth([])
# end 

# #s1 = smth([1,2,3,4,5,6])

# b1 = RemoteChannel(3)

# put!(b1, s1)

# s1 = 0; # destroy on master process. 

# r1 = remotecall(buildsm, 3)

# fetch(r1)

# remotecall_fetch(adder, 3, )

# # c1 = remotecall_fetch(adder, 3, fetch(r1), 100) # error?


# remotecall_fetch(adder, 3, s1, 200) # clearly this isn't doing what I thought.

# remotecall_fetch(adder, 3, smth([1,2,3,4,5,6]), 200)

# # this does work, but does not obviously terminate when it should. :

# remotecall_fetch(ExactVal, 3, CounterObjects(5), [11], patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0), patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0))

# remotecall_fetch(ExactVal, 4, CounterObjects(5), [1], patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0), patientcount(0.0,0.0,0.0,0.0,0.0,0.0,0.0))



# take!(b1)

# I would like to PUT objects on the subsidiary processes and refer to them by name on the main process.  




"""
`ExitComparison`
Set every firm in every zip to exit  then figure out average distances traveled.
Also prints a population and a zip code weighted average distance traveled.  

dyn = CounterObjects(1);
conf2 = [4530190 4916068 4916029 4536048 4530200 4536337 4530170 4536338 4536253];
medcounts_ex = Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }();
privcounts_ex = Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}}();

ExitComparison(conf2, medcounts_ex, privcounts_ex, 4530190)

v3 = TakeAverage(dyn, medcounts_ex, privcounts_ex, 4530190)
"""
function ExitComparison(conf::Array{Int64}, # this does take a list of fids  
                        medcounts::Dict{NTuple{9,Int64}, Dict{Int64,Array{DR,1} } }, 
                        privcounts::Dict{NTuple{9,Int64}, Dict{Int64, Array{DR,1}}},
                        special::Int64)
  d = CounterObjects(1); # call within function since this one won't get used again. 
  loc = Finder(d, special)
  chunk = [loc] # special fid. 
  k = special;  
  all_locs::Dict{Int64,Int64} = Dict{Int64,Int64}()
  neighbors::Array{Int64,1} = Array{Int64,1}()
  FindComps(d, neighbors, d.all[chunk[1]])
  push!(neighbors, chunk[1])
  CompsDict(neighbors, d, all_locs)
  state = (0,0,0,0,0,0,0,0,0)
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
  td = 0.0 # total distance*total people 
  tp = 0.0 # total people 
  tz = size(d.all[all_locs[k]].mk.m,1) # total zips 
  vtd = 0.0 # just distances - compute mean distance to a zip by dividing by tz 
  for i = 1:size(d.all[all_locs[k]].mk.m,1)
    sm = d.all[all_locs[k]].mk.m[i].pcounts.u385 + d.all[all_locs[k]].mk.m[i].pcounts.u386 + d.all[all_locs[k]].mk.m[i].pcounts.u387 + d.all[all_locs[k]].mk.m[i].pcounts.u388 + d.all[all_locs[k]].mk.m[i].pcounts.u389 + d.all[all_locs[k]].mk.m[i].pcounts.u390 + d.all[all_locs[k]].mk.m[i].pcounts.u391
    smm = d.all[all_locs[k]].mk.m[i].mcounts.u385 + d.all[all_locs[k]].mk.m[i].mcounts.u386 + d.all[all_locs[k]].mk.m[i].mcounts.u387 + d.all[all_locs[k]].mk.m[i].mcounts.u388 + d.all[all_locs[k]].mk.m[i].mcounts.u389 + d.all[all_locs[k]].mk.m[i].mcounts.u390 + d.all[all_locs[k]].mk.m[i].mcounts.u391
    ds1 = distance(d.all[all_locs[k]].mk.m[i].lat, d.all[all_locs[k]].mk.m[i].long, d.all[loc].lat, d.all[loc].long)
    td += ds1*(sm+smm)
    tp += (sm+smm)
    vtd += ds1 
  end 
  println("pop weighted avg ", td/tp, "  unweighted ", vtd/tz)
  # test...
  CleanMktDemand(d1)
  CleanMktDemand(d2)
  for k1 = 1:size(d.all[chunk[1]].mk.m,1)
    for k2 = 1:size(d.all[chunk[1]].mk.m[k1].putils, 2)
      if !(d.all[chunk[1]].mk.m[k1].putils[1,k2] == special)
        d.all[chunk[1]].mk.m[k1].putils[2,k2] = -500.0 # set ALL utilities for ALL unavailable options equal to -50.0  
      end 
    end 
    for k2 = 1:size(d.all[chunk[1]].mk.m[k1].mutils, 2)
      if !(d.all[chunk[1]].mk.m[k1].mutils[1,k2] == special)
        d.all[chunk[1]].mk.m[k1].mutils[2,k2] = -500.0 # set ALL utilities for ALL unavailable options equal to -50.0  
      end 
    end 
  end 
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
end 











function FindAllPatients(pc::patientcollection, f::Int64)
    sm::Int64 = 0
    for k1 in keys(pc.zips)
        for k2 in keys(pc.zips[k1].facilities)
            if f == k2 
                sm += sum(pc.zips[k1].mpatients)
                sm += sum(pc.zips[k1].ppatients)
            end 
        end 
    end 
    return sm 
end 








"""
`WTPchange(d::DynState)`
Take each firm.  Update the level.  Update the WTP.  Compute the change.
See also `WTPCheck` below. 

dyn = CounterObjects(1);
WTPchange(dyn)

"""
function WTPchange(d::DynState)
    outp::Array{Float64,2} = zeros(286, 7)
    for i1 in 1:size(d.all,1)
        for el in d.all[i1].mk.m # fix ALL WTP's for ALL firms 
            WTPNew(el.putils, el.pwtp)
        end  
    end
    for ix in 1:size(d.all,1)
        outp[ix,1] = d.all[ix].fid
        outp[ix,2] += d.all[ix].level
        outp[ix,3] = FindWTP(d.all[ix])
        if d.all[ix].level == 1
            outp[ix,4] = 3                                     # Checking increase given transition TO 3 FROM 1
            for zp in d.all[ix].mk.m                              
                UtilUp(zp, d.all[ix].fid, d.all[ix].level, 3)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,5] = FindWTP(d.all[ix]) 
            for zp in d.all[ix].mk.m                           # Resetting WTP                      
                UtilUp(zp, d.all[ix].fid, 3, d.all[ix].level)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,6] = 2                                     # checking increase given transition TO 2 FROM 1
            for zp in d.all[ix].mk.m                              
                UtilUp(zp, d.all[ix].fid, d.all[ix].level, 2)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,7] = FindWTP(d.all[ix]) 
            for zp in d.all[ix].mk.m                           # Resetting WTP again.                      
                UtilUp(zp, d.all[ix].fid, 2, d.all[ix].level)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
        elseif d.all[ix].level == 2
            outp[ix,4] = 3                                     # Checking increase given transition TO 3 FROM 2
            for zp in d.all[ix].mk.m                              
                UtilUp(zp, d.all[ix].fid, d.all[ix].level, 3)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,5] = FindWTP(d.all[ix]) 
            for zp in d.all[ix].mk.m                           # Resetting WTP                      
                UtilUp(zp, d.all[ix].fid, 3, d.all[ix].level)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,6] = 1                                     # checking increase given transition TO 1 FROM 2
            for zp in d.all[ix].mk.m                              
                UtilUp(zp, d.all[ix].fid, d.all[ix].level, 1)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,7] = FindWTP(d.all[ix]) 
            for zp in d.all[ix].mk.m                           # Resetting WTP again.                      
                UtilUp(zp, d.all[ix].fid, 1, d.all[ix].level)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
        elseif d.all[ix].level == 3
            outp[ix,4] = 1                                     # Checking increase given transition TO 1 FROM 3
            for zp in d.all[ix].mk.m                              
                UtilUp(zp, d.all[ix].fid, d.all[ix].level, 1)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,5] = FindWTP(d.all[ix]) 
            for zp in d.all[ix].mk.m                           # Resetting WTP                      
                UtilUp(zp, d.all[ix].fid, 1, d.all[ix].level)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,6] = 2                                     # checking increase given transition TO 2 FROM 3
            for zp in d.all[ix].mk.m                              
                UtilUp(zp, d.all[ix].fid, d.all[ix].level, 2)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end 
            outp[ix,7] = FindWTP(d.all[ix]) 
            for zp in d.all[ix].mk.m                           # Resetting WTP again.                      
                UtilUp(zp, d.all[ix].fid, 2, d.all[ix].level)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
                WTPNew(zp.putils, zp.pwtp)
            end         
        end 
    end 
    return sortrows(outp, by=x->x[2])
end 




"""
`WTPCheck(x::Array{Float64,2})`

Takes the output from `WTPChange` and computes absolute and relative changes in WTP given transfers between
levels.  Computes these changes and then prints them.  Returns nothing.   

dyn = CounterObjects(1);
wtpv = WTPchange(dyn);
WTPCheck(wtpv)

"""
function WTPCheck(x::Array{Float64,2})
    # structure of x is 
    # [fid, actual level, initial WTP, TO level, new WTP, TO level, new WTP]
    # These vectors track absolute and relative WTP differences: o(from)(to) 
    # Absolute Differences: 
    o13::Array{Float64,1} = Array{Float64,1}() # eg., from 1 to 3
    o12::Array{Float64,1} = Array{Float64,1}()
    o21::Array{Float64,1} = Array{Float64,1}()
    o23::Array{Float64,1} = Array{Float64,1}()
    o31::Array{Float64,1} = Array{Float64,1}()
    o32::Array{Float64,1} = Array{Float64,1}()
    # Relative Differences: 
    ro13::Array{Float64,1} = Array{Float64,1}() # eg., from 1 to 3
    ro12::Array{Float64,1} = Array{Float64,1}()
    ro21::Array{Float64,1} = Array{Float64,1}()
    ro23::Array{Float64,1} = Array{Float64,1}()
    ro31::Array{Float64,1} = Array{Float64,1}()
    ro32::Array{Float64,1} = Array{Float64,1}()
    # Counts 
    n1::Int64 = 0
    n2::Int64 = 0
    n3::Int64 = 0
    for i = 1:size(x, 1) # rows 
        if x[i,2] == 1
            n1 += 1
            push!(o13,  x[i,5] - x[i,3])                          # WTP at 3 - WTP at 1: should be positive.
            push!(o12,  x[i,7] - x[i,3])                          # WTP at 2 - WTP at 1: should be positive.
            push!(ro13, (((x[i,5]-x[i,3])/x[i,3]))*100 )
            push!(ro12, (((x[i,7]-x[i,3])/x[i,3]))*100 )
        elseif x[i,2] == 2
            n2 += 1
            push!(o23, x[i,5] - x[i,3])                           # this is from 2 to 3 - positive
            push!(o21, x[i,3] - x[i,7])                           # this is from 2 TO 1 - should be negative
            push!(ro23, (((x[i,5] - x[i,3])/x[i,3]))*100)   # from 2 to 3 - increase
            push!(ro21, ((x[i,3] - x[i,7])/x[i,3])*100)           # from 2 to 1 - decrease
        elseif x[i,2] == 3
            n3 += 1
            push!(o31 , x[i,3] - x[i, 5])                         # this should be negative 
            push!(o32 , x[i,3] - x[i, 7])                         # this should be negative
            push!(ro31 , ((x[i,3] - x[i, 5])/x[i,3])*100)         # this should be negative 
            push!(ro32 , ((x[i,3] - x[i, 7])/x[i,3])*100)         # this should be negative            
        else 
            println("error at ", i)
        end 
    end 
    println("Mean Changes - Absolute terms : ")
    println(" Number at Level 1: ", n1)
    println("Changing 1 to 3 - Absolute Change in WTP: ", round(mean(o13), 3), " Standard Deviation: ", round(std(o13), 3), " Extrema: ", extrema(o13), " " )
    println("Changing 1 to 2 - Absolute Change in WTP: ", round(mean(o12), 3), " Standard Deviation: ", round(std(o12), 3), " Extrema: ", extrema(o12), " " )
    println(" Number at Level 2: ", n2)
    println("Changing 2 to 1 - Absolute Change in WTP: ", round(mean(o21), 3), " Standard Deviation: ", round(std(o21), 3), " Extrema: ", extrema(o21), " " )
    println("Changing 2 to 3 - Absolute Change in WTP: ", round(mean(o23), 3), " Standard Deviation: ", round(std(o23), 3), " Extrema: ", extrema(o23), " " )
    println(" Number at Level 3: ", n3)
    println("Changing 3 to 1 - Absolute Change in WTP: ", round(mean(o31), 3), " Standard Deviation: ", round(std(o31), 3), " Extrema: ", extrema(o31), " " )
    println("Changing 3 to 2 - Absolute Change in WTP: ", round(mean(o32), 3), " Standard Deviation: ", round(std(o32), 3), " Extrema: ", extrema(o32), " " )
    println("Mean Changes - Relative Terms :")
    println(" Number at Level 1: ", n1)
    println("Changing 1 to 3 - Relative Change in WTP: +", round(mean(ro13), 2), "% Standard Deviation: ",round(std(ro13), 3), " Extrema: ", extrema(ro13), " " )
    println("Changing 1 to 2 - Relative Change in WTP: +", round(mean(ro12), 2), "% Standard Deviation: ",round(std(ro12), 3), " Extrema: ", extrema(ro12), " " )
    println(" Number at Level 2: ", n2)
    println("Changing 2 to 1 - Relative Change in WTP: -", round(mean(ro21), 2), "% Standard Deviation: ",round(std(ro21), 3), " Extrema: ", extrema(ro21), " " )
    println("Changing 2 to 3 - Relative Change in WTP: +", round(mean(ro23), 2), "% Standard Deviation: ",round(std(ro23), 3), " Extrema: ", extrema(ro23), " " )
    println(" Number at Level 3: ", n3)
    println("Changing 3 to 1 - Relative Change in WTP: -", round(mean(ro31), 2), "% Standard Deviation: ",round(std(ro31), 3), " Extrema: ", extrema(ro31), " " )
    println("Changing 3 to 2 - Relative Change in WTP: -", round(mean(ro32), 2), "% Standard Deviation: ",round(std(ro32), 3), " Extrema: ", extrema(ro32), " " )
    return nothing 
end 


"""
`WTPchangeALL()`
Compute all possible transitions between states. 

WTPchangeALL(dyn)

"""
function WTPchangeALL(d::DynState)
    outp::Array{Any,2} = Array{Any,2}(286, 15)
    for i1 in 1:size(d.all,1)
        for el in d.all[i1].mk.m # fix ALL WTP's for ALL firms 
            WTPNew(el.putils, el.pwtp)
        end  
    end
    for ix in 1:size(d.all,1)
        outp[ix,1] = d.all[ix].fid
        outp[ix,2] = d.all[ix].level
        outp[ix,3] = FindWTP(d.all[ix])
        # Set all firms down to WTP at 1.  
        for zp in d.all[ix].mk.m                              
            UtilUp(zp, d.all[ix].fid, d.all[ix].level, 1)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
            WTPNew(zp.putils, zp.pwtp)
        end 
        outp[ix,4] = (1,3)                                  # Checking increase given transition TO 3 FROM 1
        for zp in d.all[ix].mk.m                              
            UtilUp(zp, d.all[ix].fid, 1, 3)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
            WTPNew(zp.putils, zp.pwtp)
        end 
        outp[ix,5] = FindWTP(d.all[ix]) 
        # LEAVE utility at level 3, then set down to 2.
        outp[ix,6] = (3,2)                                     # checking downgrade 3 to 2
        for zp in d.all[ix].mk.m                              
            UtilUp(zp, d.all[ix].fid, 3, 2)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
            WTPNew(zp.putils, zp.pwtp)
        end 
        outp[ix,7] = FindWTP(d.all[ix]) 
        # LEAVE utility at 2, the set up to 3
        outp[ix,8] = (2,3)                                     # Checking increase given transition TO 3 FROM 2
        for zp in d.all[ix].mk.m                              
            UtilUp(zp, d.all[ix].fid, 2, 3)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
            WTPNew(zp.putils, zp.pwtp)
        end 
        outp[ix,9] = FindWTP(d.all[ix]) 
        # LEAVE utility at 3, then set to 1
        outp[ix,10] = (3,1)                                     # checking increase given transition TO 3 FROM 1
        for zp in d.all[ix].mk.m                              
            UtilUp(zp, d.all[ix].fid, 3, 1)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
            WTPNew(zp.putils, zp.pwtp)
        end 
        outp[ix,11] = FindWTP(d.all[ix]) 
        # LEAVE utility at 1
        outp[ix,12] = (1,2)                                     # Checking increase given transition TO 1 FROM 2
        for zp in d.all[ix].mk.m                              
            UtilUp(zp, d.all[ix].fid, 1, 2)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
            WTPNew(zp.putils, zp.pwtp)
        end 
        outp[ix,13] = FindWTP(d.all[ix]) 
        # LEAVE utility at 2
        outp[ix,14] = (2,1)                                     # checking increase given transition TO 2 FROM 3
        for zp in d.all[ix].mk.m                              
            UtilUp(zp, d.all[ix].fid, 2, 1)  # UtilUp(c::cpats, fid::Int64, actual::Int64, current::Int64)
            WTPNew(zp.putils, zp.pwtp)
        end 
        outp[ix,15] = FindWTP(d.all[ix])        
    end 
    return sortrows(outp, by=x->x[2])
end 









"""
`DUtilCheck`
See whether the utilities are updated correctly for "shortrecs"
"""
function DUtilCheck(h::simh)
    println("Facilities")
    for el in h.mk.m 
        println(el.zp, " ", size(unique(el.putils[1,:]))[1], "  ", size(unique(el.mutils[1,:]))[1])
    end 


    println("Own Utility Update Check: ")
    for el in h.mk.m 
        println(el.zp, "  ", el.putils[2, findin(el.putils[1,:], h.fid)])
    end 
    h.level = 3
    h.tbu = true 
    UpdateDUtil(h)
    println("********")
    for el in h.mk.m 
        println(el.zp, "  ", el.putils[2, findin(el.putils[1,:], h.fid)])
    end
    h.level = 1 
    h.tbu = true 
    UpdateDUtil(h) 
end 




""" 
`CompareResults(d::DynState, m::Array{Any,2})`

Finds the actual configuration of the market in f and then
returns the values in the array m (containing sim results)
which correspond to that state.  

dynresults2017-09-09-10-53-56.csv
dynresults2017-09-15-07-00-16.csv

dres1 = readcsv("dynresults2017-09-09-10-53-56.csv", header = false)

dres2 = readdlm("dynresults2017-09-15-07-00-16.csv", header = false)

"""
function CompareResults(d::DynState, m::Array{Any,2}, f::Int64)
    ix = 0
    col =  # indicates column in which state tuples are located
    vals = # indicates column in which values located.
    for i = 1:size(d.all,1)
        if d.all[i].fid == f 
            ix = i # returns the index of i 
        end     
    end 
    cf1, cf2, cf3 = CNStoTupleLevel(dyn.all[ix].cns)
    for i = 1:size(m,1)
        if m[i,1] == f 
            if m[i,col] == cf1 || m[i,col] == cf2 || m[i,col] == cf3
                println("row, col ", i," ",col, "  ", m[i,col])
                println(m[i,vals])
            end 
        end 
    end 
end 










#=
This computation relates to determining the market sizes of individual hositals by checking the total number of births in 
all zip codes attached to that hospital every year.  It is merged with some data from Stata.  Specifically, see TX Patient Uncertainty.do

# =#
# mx = 0;
# for el in dyn.all 
#     if size(el.mk.m,1) > mx 
#         mx = size(el.mk.m,1)
#     end 
# end 
# println("Max zips = ", mx)

# zps = zeros(Int64, mx+1, size(dyn.all,1) ); # NB - each COLUMN is an fid 

# for el in 1:size(dyn.all,1) # COLUMNS
#     zps[1, el] = dyn.all[el].fid  
#     for z in 1:size(dyn.all[el].mk.m,1) # ROWS 
#         zps[z+1,el] = dyn.all[el].mk.m[z].zp # zip 
#     end 
# end 
# zps = transpose(zps)
# writecsv( "/Users/austinbean/Desktop/dynhosp/hospzips.csv", zps)




# #=

# typetest 


# =#

# mutable struct WT 
#     w::Array{Float64,1}
# end 

# struct WTPee
#     w::Array{Float64, 1}
# end

# w1 = WTPee(Array{Float64,1}(50))
# w2 = WT(Array{Float64,1}(50))



# for i = 1:10
# @time for k1 = 1:50 
#     w1.w[k1] = 3.0
# end 
# for k = 1:50 
#     w1.w[k] = 0.0
# end 
# end 

# for i = 1:10
# @time for k1 = 1:50 
#     w2.w[k1] = 3.0
# end 
# for k = 1:50 
#     w2.w[k] = 0.0
# end 
# end 



# struct w1
#     d::Dict{Int64,Int64}
# end 

# ab1 = WTPee(Array{Float64,1}())

# ab2 = WTPee(Array{Float64,1}(50))

# @time for el =1:50
#     push!(ab1.w, rand())
# end 

# @time for el =1:50 
#     ab2.w[el] = rand()
# end 

# for j = 1:10
# @time for i = 1:50
#  ab2.w[i] = i 
# end 
# for k = 1:50
#     ab2.w[k] = 0.0
# end 
# end 


# ab = w1(Dict{Int64,Int64}())


# for el in Texas.ms
#     for k in el.config 
#         levl = (-1,-1)
#          if k.level == 1
#            levl = (0,0)
#          elseif k.level == 2
#            levl = (1,0)
#          elseif k.level == 3
#            levl = (0,1)
#          end
#         n1, n2, n3 = MktSize(k.neigh)
#         nv = [k.neigh.level105; k.neigh.level205; k.neigh.level305; k.neigh.level1515; k.neigh.level2515; k.neigh.level3515; k.neigh.level11525; k.neigh.level21525; k.neigh.level31525 ]
#         try a1 = logitest(levl, n1, n2, n3, nv)
#         catch er1
#             println(er1, "the error")
#             if isa(er1,ProjectModule.ValueException)
#                 println("didn't work:")
#                 println(levl)
#                 println(n1, " ", n2, " ", n3)
#                 println(nv)
#                 println(typeof(nv))
#             else 
#                 # do nothing
#             end 
#         end 
#     end 
# end 



# #=
# import Base.+
# Base.+(a1::String, a2::String) = string(a1,a2)
# =#





# # TODO - must get the output of ResultsOut and ResultsOutVariant into the @parallel.  

# Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);

# outp1, outp2 = @sync @parallel (+) for j = 1:MCcount
#     # probably need a new function to make this work.  

# end 

# ut = [0.1, 0.2, 0.3, 0.4, 0.5];
# fi = [111, 222, 333, 444, 19];
# ta = [0.0, 0.0, 0.0, 0.0, 0.0];
# function UMap(utils::Array{Float64,1},
#               fids::Array{Int64,1},
#               temparr::Array{Float64,1})::Int64
#   const dist_μ::Int64 = 0
#   const dist_σ::Int64 = 1
#   const dist_ξ::Int64 = 0
#   d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
#   return fids[indmax(utils+rand!(d, temparr))]
# end



# function NUMap(utils::Array{Float64,1},fids::Array{Int64,1})::Int64
#   const dist_μ::Int64 = 0
#   const dist_σ::Int64 = 1
#   const dist_ξ::Int64 = 0
#   d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
#   maxm::Int64 = 1                   # this is an index to a fid.  Will record max index.
#   maxu::Float64 = 0.0               # holds the value of max utility 
#   tem::Float64 = 0.0                # holds interim utility value
#   for i = 1:size(utils,1)
#     tem = utils[i]+rand(d)
#     if tem>maxu  
#         maxm = i                    # replace index if greater
#         maxu = tem                  # replace max util 
#     end 
#   end 
#   return fids[maxm]::Int64          # return fid of max util value.  
# end 

# UMap(ut,fi,ta)
# NUMap(ut,fi)
# srand(1)
# for i = 1:10
#     println(UMap(ut, fi, ta))
# end 
# srand(1)
# for i = 1:10
#     println(NUMap(ut,fi))
# end 

# # Some Benchmarks related to PSim()
# # this is a pre-change in UMap benchmark.
# @benchmark GenPChoices($patients, $d1, $arry1)
# BenchmarkTools.Trial:
#   memory estimate:  33.44 MiB
#   allocs estimate:  644479
#   --------------
#   minimum time:     143.279 ms (1.80% GC)
#   median time:      147.961 ms (3.44% GC)
#   mean time:        148.439 ms (2.72% GC)
#   maximum time:     153.695 ms (3.65% GC)
#   --------------
#   samples:          34
#   evals/sample:     1

# # 
# @benchmark ChoiceVector($patients.zips[78759].pdetutils, $dic1, $inpt, $patients.zips[78759].ppatients)
# BenchmarkTools.Trial:
#   memory estimate:  19.20 KiB
#   allocs estimate:  1156
#   --------------
#   minimum time:     298.219 μs (0.00% GC)
#   median time:      301.069 μs (0.00% GC)
#   mean time:        310.506 μs (1.49% GC)
#   maximum time:     6.268 ms (94.06% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1


# # this is a benchmark for the OLD UMap.  
#     @benchmark UMap($ut, $fi, $ta)
#         BenchmarkTools.Trial:
#           memory estimate:  128 bytes
#           allocs estimate:  1
#           --------------
#           minimum time:     324.303 ns (0.00% GC)
#           median time:      327.463 ns (0.00% GC)
#           mean time:        346.376 ns (1.78% GC)
#           maximum time:     6.011 μs (83.58% GC)
#           --------------
#           samples:          10000
#           evals/sample:     231

# # this will be the new UMap - 10% speedup.
# @benchmark NUMap($ut, $fi)
# BenchmarkTools.Trial:
#   memory estimate:  0 bytes
#   allocs estimate:  0
#   --------------
#   minimum time:     283.248 ns (0.00% GC)
#   median time:      288.248 ns (0.00% GC)
#   mean time:        297.208 ns (0.00% GC)
#   maximum time:     3.209 μs (0.00% GC)
#   --------------
#   samples:          10000
#   evals/sample:     298

# @benchmark DV($patients.zips[78702].pdetutils)
# BenchmarkTools.Trial:
#   memory estimate:  400 bytes
#   allocs estimate:  4
#   --------------
#   minimum time:     349.616 ns (0.00% GC)
#   median time:      357.720 ns (0.00% GC)
#   mean time:        388.092 ns (6.41% GC)
#   maximum time:     8.801 μs (92.44% GC)
#   --------------
#   samples:          10000
#   evals/sample:     211




# @benchmark WriteWTP($WTPMap(patients, Texas), $Texas, 1)
# BenchmarkTools.Trial:
#   memory estimate:  2.56 MiB
#   allocs estimate:  14741
#   --------------
#   minimum time:     6.164 ms (0.00% GC)
#   median time:      7.412 ms (0.00% GC)
#   mean time:        7.681 ms (3.55% GC)
#   maximum time:     12.744 ms (31.44% GC)
#   --------------
#   samples:          650
#   evals/sample:     1

# @benchmark WTPMap($patients, $Texas)
# BenchmarkTools.Trial:
#   memory estimate:  2.50 MiB
#   allocs estimate:  10737
#   --------------
#   minimum time:     5.307 ms (0.00% GC)
#   median time:      5.828 ms (0.00% GC)
#   mean time:        6.185 ms (3.11% GC)
#   maximum time:     10.855 ms (26.12% GC)
#   --------------
#   samples:          807
#   evals/sample:     1

# @benchmark NewHospDict($Texas)
# BenchmarkTools.Trial:
#   memory estimate:  42.53 KiB
#   allocs estimate:  304
#   --------------
#   minimum time:     18.695 μs (0.00% GC)
#   median time:      21.220 μs (0.00% GC)
#   mean time:        36.018 μs (14.99% GC)
#   maximum time:     5.607 ms (97.76% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1

# @benchmark EntryProcess($Texas.mkts[48453], 1, 50)
# BenchmarkTools.Trial:
#   memory estimate:  441 bytes
#   allocs estimate:  3
#   --------------
#   minimum time:     20.789 μs (0.00% GC)
#   median time:      76.908 μs (0.00% GC)
#   mean time:        79.714 μs (0.17% GC)
#   maximum time:     150.223 μs (0.00% GC)
#   --------------
#   samples:          74
#   evals/sample:     858



# # getting rid of the patient counts:  

# Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
# patients = CreateZips(ProjectModule.alldists, Texas);
# FillPatients(patients , ProjectModule.pinsured, ProjectModule.pmedicaid)

# "/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 Medicaid Individual Choices.csv"



#=

=#