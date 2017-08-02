

"""
`NonZeros(d::Dict{Int64,Dict{NTuple{10,Int64},Float64}})
"""
function NonZeros(d::Dict{Int64,Dict{NTuple{10,Int64},Float64}})
  for k1 in keys(d)
    println("For ", k1)
    ct::Int64 = 0
    for k2 in keys(d[k1])
      if d[k1][k2]>0
        ct += 1
        println(k2, " ", d[k1][k2])
      end 
    end 
    println("total nonzero: ", ct)
  end 
end 


"""
`MktSize(d::DynState)`
How many patients of each type are there in each firm's market?
"""
function MktSize(d::DynState)
  emp::Array{Int64,1} = Array{Int64,1}()
  for el in d.all 
    mp::Int64 = 0
    pp::Int64 = 0
    for z in el.mk.m  
      pp += sum(z.pcounts)
      if sum(z.pcounts) == 0
        push!(emp,z.zp)
      end 
      mp += sum(z.mcounts)
      if sum(z.mcounts) == 0
        push!(emp, z.zp)
      end 
    end 
    println(el.fid, " private: ", pp, " medicaid: ", mp)
  end 
  println("Empty zips: ", size(unique(emp)))
  for z in unique(emp)
    println(z)
  end 
end 




"""
`MarketPrint(mkt::Market)`
Prints the elements of the market record: name, neighbors, choice probabilities.
For debugging purposes to make sure things look right.
"""
function MarketPrint(mkt::Market)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  println(mkt.fipscode)
  for el in mkt.config
    println("*******************")
    println(el.name)
    println(el.neigh)
    println(el.hood)
    println(el.chprobability)
  end
end





"""
`NeighborsPrint(mkt::Market)`
Prints the name and neighbors of every hospital in the market `mkt`
"""
function NeighborsPrint(mkt::Market)
  for el in mkt.config
    println(el.name, " ", el.neigh)
  end
end




"""
`FacPrint(hosp::hospital)`
Prints out a hospital facility record: name, fid, fips, level, neighbors.
"""
function FacPrint(hosp::hospital)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  println(hosp.name)
  println(hosp.fid)
  println("Fips: ", hosp.fipscode)
  println("Level: ", hosp.level)
  println("Neighbors: ", hosp.hood)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
end




"""
`PrintZip(zi::zip)`
Prints the fid and the name of the facilities attached to the zips.
"""
function PrintZip(zi::zipcode)
  for el in keys(zi.facilities)
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
    println(el, "  ", zi.facilities[el].name)
  end
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
    println("385 ", zi.ppatients.count385, " ", zi.mpatients.count385)
    println("386 ", zi.ppatients.count386, " ", zi.mpatients.count386)
    println("387 ", zi.ppatients.count387, " ", zi.mpatients.count387)
    println("388 ", zi.ppatients.count388, " ", zi.mpatients.count388)
    println("389 ", zi.ppatients.count389, " ", zi.mpatients.count389)
    println("390 ", zi.ppatients.count390, " ", zi.mpatients.count390)
    println("391 ", zi.ppatients.count391, " ", zi.mpatients.count391)
end



"""
`ZeroFind(mat::Array{Float64,2})`
looks for hospitals which were never demanded - there are some which have rows of zeros
will want to pop those out of the results matrix.
"""
function ZeroFind(mat::Array{Float64,2})
  zers = zeros(size(mat,1)) # vector of zeros for output
  for i = 1:size(mat,1)
    zercnt = 0
    for j = 1:size(mat, 2) # don't want to count all of them - fix this.
      if mat[i,j] == 0.0
        zercnt += 1
      end
    end
    zers[i] = zercnt
  end
  return zers
end


"""
`CheckPats(pats::patientcollection)`
This will compute the total sum of all patients in all zips in the patientcollection.
This is just to check that all are being added as expected.
"""
function CheckPats(pats::patientcollection)
  count::Int64 = 0
  for k in keys(pats.zips)
    count += sum(pats.zips[k].mpatients)
    count += sum(pats.zips[k].ppatients)
  end
  return count
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
`ExitUpdate(h::simh, l::Int64)::Array{Float64,2}`
Checks that the update of h to level l works:
- Records utils across zips
- Changes level to l, updates utils, records those 
- Changes level back, resets utils, records them 
- Return results by zip in an array.  

- Further test... change neighbors too.   UGH.  

"""
function ExitUpdate(h::simh, l::Int64)::Array{Float64,2}
  dim1 = size(h.mk.m,1)
  original = dyn.all[1].level
  outp::Array{Float64,2}=Array{Float64,2}(dim1, 4)
  for el in 1:size(h.mk.m,1)
     outp[el,1] = h.mk.m[el].zp
     for f1 in 1:size(h.mk.m[el].putils,2)
        if h.fid == h.mk.m[el].putils[1,f1]
           outp[el,2] = h.mk.m[el].putils[2,f1] 
        end 
     end 
  end 
  # do the update 
  h.actual = h.level 
  h.level = l
  UpdateD(h)
  # check the results.
  for el in 1:size(h.mk.m,1)
    outp[el,1] = h.mk.m[el].zp
    for f1 in 1:size(h.mk.m[el].putils,2)
      if h.fid == h.mk.m[el].putils[1,f1]
        outp[el,3] = h.mk.m[el].putils[2,f1] 
      end 
    end 
  end 
  # return to original
  println("original: ", original)
  h.actual = h.level 
  h.level = original 
  UpdateD(h)
  for el in 1:size(h.mk.m,1)
    outp[el,1] = h.mk.m[el].zp
    for f1 in 1:size(h.mk.m[el].putils,2)
      if h.fid == h.mk.m[el].putils[1,f1]
        outp[el,4] = h.mk.m[el].putils[2,f1] 
      end 
    end 
  end
  return outp
end 

"""
`UtilCCheck`
Check whether the U up and down functions correctly restore to the original state.

dyn = CounterObjects(5);
dyn2 = CounterObjects(5);
UtilCCheck(dyn, dyn2)


this should maybe call a new state, map it out, do the update, then check all the utilities across the values.  
That can happen just once, let's say.  

"""
function UtilCCheck(d1::DynState, d2::DynState)
  # what do I want this to check?  Create a state, write out some different utilities, 
  # change them using UtilUp and UtilDown, then figure out where the fuckup is coming.
  # This is facilitated by comparing the utility to the original value in d2.
  outvals::Dict{ Int64, Dict{NTuple{10, Int64},  Float64} } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  tempvals::Dict{ Int64, Dict{NTuple{10, Int64}, Float64}  } = Dict{ Int64, Dict{NTuple{10, Int64}, Float64 } }()
  totest::Dict{Int64,Bool} = Dict{Int64,Bool}()
  el = 11 # testing this on a specific facility. 
  chunk = [11]; 
  k = 2910645; # fid of dyn.all[el]
  st_dict::Dict{Int64,NTuple{9,Int64}} = Dict{Int64,NTuple{9,Int64}}()
  neighbors::Array{Int64,1} = Array{Int64,1}()                                          
  nfds::Array{Int64,1} = Array{Int64,1}()
  all_locs::Dict{Int64,Int64} = Dict{Int64, Int64}()
  for el in chunk                                                                       # goal of this loop is to: set up the dictionaries containing values with entries for the fids.  
    FindComps(d1, neighbors, d1.all[el])                                                  # these are addresses of competitors in D.all 
    NFids(d1, nfds, d1.all[el])                                                           # records the fids of neighbors only, as fids.  
    push!(neighbors, el)                                                                # add the location of the firm in chunk
    if !haskey(outvals, d1.all[el].fid)                                                  # add an empty dict IF there isn't already an entry.
      outvals[d1.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    end 
    tempvals[d1.all[el].fid] = Dict{NTuple{10, Int64}, Float64 }()
    CompsDict(neighbors, d1, all_locs)                                                   # now this is Dict{Fid, Location}
    StateRecord(all_locs, d1, st_dict)                                                   # returns the restricted state.
    for el2 in neighbors                                                                # adds keys for the neighbors to the temp dict. 
      outvals[d1.all[el2].fid] = Dict{NTuple{10, Int64}, Float64}()
      tempvals[d1.all[el2].fid] = Dict{NTuple{10, Int64},Float64}() 
      totest[d1.all[el2].fid] = false                                                    # initialized to FALSE - not converged. 
      StateEnumerate(TupletoCNS(st_dict[d1.all[el2].fid]), outvals[d1.all[el2].fid]) 
      StateEnumerate(TupletoCNS(st_dict[d1.all[el2].fid]), tempvals[d1.all[el2].fid])
    end 
    if !haskey(outvals, d1.all[el].fid)
      StateEnumerate(d1.all[el].cns, outvals[d1.all[el].fid])
      StateEnumerate(d1.all[el].cns, tempvals[d1.all[el].fid])                            # this does NOT need starting values.  
    end 
    if haskey(outvals, d1.all[el].fid)
      DictCopy(tempvals, outvals, 1.0)                                                  # if there is an entry for the value, copy FROM outvals TO tempvals.  
    end                                
    totest[d1.all[el].fid] = false                                                       # all facilities to do initially set to false.  
  end
  #=
  Other States to test:
  Tuple{Int64,Int64}[(3396189, 3), (3396327, 3), (3390720, 1), (3396057, 2)]
  Tuple{Int64,Int64}[(3396189, 1), (3396327, 1), (3390720, 2), (3396057, 2)]
  Tuple{Int64,Int64}[(3396189, 1), (3396327, 1), (3390720, 1), (3396057, 2)]

  [(3396189, 1) (3396327, 1) (3390720, 1) (3396057, 1); (3396189, 2) (3396327, 1) (3390720, 1) (3396057, 1); (3396189, 3) (3396327, 1) (3390720, 1) (3396057, 1); (3396189, 999) (3396327, 1) (3390720, 1) (3396057, 1); (3396189, 1) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 2) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 3) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 999) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 1) (3396327, 3) (3390720, 1) (3396057, 1); (3396189, 2) (3396327, 3) (3390720, 1) (3396057, 1)]
  =#
  altstates =   [(3396189, 1) (3396327, 999) (3390720, 1) (3396057, 1); (3396189, 2) (3396327, 999) (3390720, 999) (3396057, 1); (3396189, 3) (3396327, 1) (3390720, 1) (3396057, 1); (3396189, 999) (3396327, 1) (3390720, 1) (3396057, 1); (3396189, 1) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 2) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 3) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 999) (3396327, 2) (3390720, 1) (3396057, 1); (3396189, 1) (3396327, 3) (3390720, 1) (3396057, 1); (3396189, 2) (3396327, 3) (3390720, 1) (3396057, 1)]
  for k in keys(totest)                                                              
    if !totest[k]  
      for r in 1:size(altstates,1)
        #println(altstates[r,:])
        st_dict[k] = GiveState(d1, chunk, all_locs, altstates[r,:], d1.all[el].cns)
        MapCompState(d1, all_locs, chunk, FindFids(d1, chunk), altstates[r,:]) 
        println("Pre-audit: ")
        DynAudit(d1, d2)
        ExactChoice(tempvals, outvals, all_locs, st_dict, k, all_locs[k], p1, p2, d1; messages = false) 
        ResetCompState(d1, all_locs, chunk, FindFids(d1, chunk), altstates[r,:])
        println("Post audit: ")
        DynAudit(d1, d2)
      end 
    end 
  end 
end


"""
`DynAudit`
Checks the creation of dynstate, especially the utility and WTP levels.
The condition in isapprox should nearly always be triggered for all firms and zips, 
but differences should be the same as those appearing in `UtilUp`

dyn = CounterObjects(5);
DynAudit(dyn)

d.all[el].mk.m[z].putils[2,i], d2.all[el].mk.m[z].putils[2,i]
  # for h in d.all 
  #   UpdateUCheck(h)
  # end 

  dyn = CounterObjects(5);
  d2 = CounterObjects(5);
  dyn.all[1].mk.m[1].putils[2,3] += 20.0
  DynAudit(dyn, d2)

  diffs = [-0.57239721 0.57239721 1.3499395 0.77754229 0.31972186 -0.31972186 1.18599 0.866268]
  [ abs(i-j) for i in diffs, j in diffs]

 
"""
function DynAudit(d::DynState, d2::DynState)
  diffs::Array{Float64} = [-0.57239721 0.57239721 1.3499395 0.77754229 0.31972186 -0.31972186 1.18599 0.866268]
  for el in 1:size(d.all,1) 
    for z in 1:size(d.all[el].mk.m,1) 
      for i = 1:size(d.all[el].mk.m[z].putils,2)
        if !isapprox(d.all[el].mk.m[z].putils[2,i], d2.all[el].mk.m[z].putils[2,i], atol=10e-5)
          if !CollectApprox(abs(d.all[el].mk.m[z].putils[2,i]- d2.all[el].mk.m[z].putils[2,i]),diffs)
            if abs(d.all[el].mk.m[z].putils[2,i]- d2.all[el].mk.m[z].putils[2,i]) > 20 # catches exit/returns which are not fixed.  
              println(d.all[el].fid, " ", d.all[el].mk.m[z].zp, " ", convert(Int64,d.all[el].mk.m[z].putils[1,i])," ", d.all[el].level, " ", d.all[el].actual  , " ", d.all[el].mk.m[z].putils[2,i]- d2.all[el].mk.m[z].putils[2,i] )
            end 
          end 
        end 
      end 
    end 
  end 
end 


"""
`CollectApprox`
Is the element approximately equal to one element in the collection?

test1 = [0.57239721 1.3499395 0.77754229 0.31972186 1.18599 0.866268]
CollectApprox(1.185991400000000, test1)


"""
function CollectApprox(x, coll; toler::Float64 = 10e-4)
  res::Bool = false 
  for el in coll 
    res = res||isapprox(x, el; rtol = toler)
    if res 
      break 
    end 
  end 
  return res 
end 


"""
`UpdateUCheck(h::simh)`
This checks that the utility update actually worked.  Apply in ExactChoice.  
Find the WTP shares and print them too. 

TODO - there is a utility computation problem.  The lowest value is way too low.  
This is not an initialization problem.  This is a change problem.  
Clearly the problem is with the util up and down functions and the way they are called.  

"""
function UpdateUCheck(h::simh)
  #   WTPNew(inparr, temparr) # updates temparr - it should be with shares given by utilities. 
  #   ArrayZero(temparr)
  temparr = zeros(2,12) # for WTP.  
  for el in h.mk.m 
    u1, iu1 = findmax(el.putils[2,:])
    u2, iu2 = findmin(el.putils[2,:])
    WTPNew(el.putils, temparr) # what array dimensions? 
    m1, i1 = findmax(temparr[2,:]) 
    #m2, i2 = findmin(temparr[2,:]) # not informative since frequently zero. 
    for j = 1:size(el.putils,2)
      if el.putils[1,j] == h.fid 
       println(el.zp, "  ", h.fid, " UTIL: ",  round(el.putils[2,j],4), " WTP: ", round(temparr[j],4), " MAX WTP: ", round(m1,4), " FAC: ", el.putils[1,i1],  " MAX U: ", round(u1,4), " FAC: ", el.putils[1,iu1], " MIN U: ", round(u2,4), " FAC: ", el.putils[1,iu2] )
      end 
    end
    ArrayZero(temparr)
  end 
end 



"""
`CpatsUChange(c::cpats, f::Int64)`

take a fid and a cpat, find the fid in putils and mutils, print out the fid and the utility value.
"""
function CpatsUChange(c::cpats, f::Int64)
  for el in 1:size(c.putils,2)
    if c.putils[1,el] == f 
      print(f, " ", c.putils[2,el])
    end 
  end 
  for el in 1:size(c.mutils,2) 
    if c.mutils[1,el] == f 
      println(" ", c.mutils[2,el])
    end 
  end 
end


