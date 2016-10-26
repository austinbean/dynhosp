# Now solve the dynamic model itself.


# Can put these into a market and then into a state.

# use the type neighbors !

using Distributions
using DataFrames
type hitcount # this type will record visits to a state-action pair
  conf::neighbors
  visits::Dict{Int64,Int64}
end

type history
  path::Dict{neighbors, hitcount}
  totalcount::Int64 # records total number of iterations
end


# Initialize empty:  hithist = history(Dict{neighbors, hitcount}(), 0)

type hstate
  #NB: track these for the level !
  conf::neighbors
  wvalues::Array{Float64,1}
  actprobs::Array{Float64,1}
end


#NB: here include a bool for "update" - that will track that some facility should be updated among the patients.
type shortrec<:ProjectModule.Fac
  # needs to take location.  must measure distances.  need to update all of these every time, probably.
  fid::Int64
  lat::Float64
  long::Float64
  level::Int64
  truelevel::Int64
  ns::neighbors
  choices::Array{Int64, 2}
  chprobs::WeightVec
  tbu::Bool
end

#NB: here include a boolean that some facility should be updated among the patients.
type simh<:ProjectModule.Fac
  fid::Int64
  lat::Float64
  long::Float64
  level::Int64
  actual::Int64
  cns::neighbors # must know what current neighbors look like.
  visited::Dict{neighbors, hstate} #possible given "isequal" and "hash" extended for "neighbors"
  ns::Array{shortrec, 1}
  exit::Bool
  tbu::Bool
end


type DynState # this should hold a collection of ALL of the hospitals, for the approximation.
  all::Array{simh, 1}
end


#NB: It may be easiest eventually to put this collection directly in the simh.

type cpats
  zp::Int64
  lat::Float64
  long::Float64
  putils::Array{Float64,2}
  mutils::Array{Float64,2}
  pcounts::patientcount
  mcounts::patientcount
end

type cmkt
  fid::Int64
  m::Array{cpats,1}
end


"""
`DynStateCreate(Tex::EntireState)`
Create the dynamic records from the existing state, don't bother doing it from scratch.
And these don't have to be organized by zip.
Use the EntireState from the equilibrium simulation, not the first counterfactual.
Make using
TexasEq = MakeNew(ProjectModule.fips, ProjectModule.data05);
dyn = DynStateCreate(TexasEq);
"""
function DynStateCreate( Tex::EntireState )
  outp = DynState(Array{simh,1}())
  for k1 in keys(Tex.mkts)
    for hk in keys(Tex.mkts[k1].collection)
      newsimh = simh(Tex.mkts[k1].collection[hk].fid,
                     Tex.mkts[k1].collection[hk].lat,
                     Tex.mkts[k1].collection[hk].long,
                     Tex.mkts[k1].collection[hk].level,
                     Tex.mkts[k1].collection[hk].level,
                     neighbors(0,0,0,0,0,0,0,0,0),
                     Dict{neighbors,hstate}(),
                     Array{shortrec,1}(),
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
      push!(outp.all, newsimh)
    end
  end
  return outp
end


"""
`DynPatients(p::patientcollection, f::fid)`
This will take a collection of patients and create a cmkt, which is a vector of
cpats.  That is, it will take every zip for which `f` is an option, then create
the collection of patients for those zips.
"""
function DynPatients(p::patientcollection, f::Int64 )
  outp = cmkt(f, Array{cpats,1}())
  zpc = PatientFind(p, f)
  for el in zipc
    push!(outp.m, cpats(el,
                  p.zips[el].lat,
                  p.zips[el].long,
                  DetUtils(p.zips[el]), # write this - use hcat
                  DetUtils(p.zips[el]),
                  patients...,
                  patients... )
  end
  return outp
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
`Simulate(d::DynState)`
This should take a dynamic state, as generated above in `DynStateCreate` and do a simulation.
What does that mean?
- Start in some state for each firm
- Choose actions by the other firms in the simh.ns
- Compute the demand and the return for the simulating hospital
- Write out the return to this, updating the estimate of the return at that state.
"""
function Simulate(d::DynState)



end




# The "Payoff" function in counter 1 gets the profit

"""
`ComputeR(hosp::simh, ppats::Dict{Int64, ProjectModule.patientcount}, mpats::Dict{Int64, ProjectModule.patientcount}, Tex::EntireState, wtp::Dict{Int64,Float64} )`
Computes the return (current profit + expected continuation) for each hospital in the state.
"""
function ComputeR(hosp::simh,
                  ppats::Dict{Int64, ProjectModule.patientcount},
                  mpats::Dict{Int64, ProjectModule.patientcount},
                  Tex::EntireState,
                  wtp::Dict{Int64,Float64} ;
                  disc::Float64 = 0.95)
  #NB: all of this can be updated in place, I think.
  for

    Payoff() + disc*(dot(hosp.wvalues, hosp.actprobs) + eulergamma - dot(log(hosp.actprobs),hosp.actprobs))

  end
end

"""
`PolicyUpdate(hosp::simh, neww::Array{Float64,1})`
Takes an array of the new W values and maps back to the simh action probabilities.
"""
function PolicyUpdate(hosp::simh, neww::Array{Float64,1})
  hosp.actprobs = exp(neww - maximum(neww))/sum(exp(neww-maximum(neww)))
end


"""
`UpdateW(hosp::simh, ret::Float64, stcounter::Int64)`
"""
function UpdateW(hosp::simh, ret::Float64, stcounter::Int64)


end



function CheckConvergence()


end
