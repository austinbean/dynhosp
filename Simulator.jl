using DataFrames
using DataArrays
using Distributions

include("/Users/austinbean/Desktop/dynhosp/LogitEst.jl")
include("/Users/austinbean/Desktop/dynhosp/Distance.jl")

T = 100;
neighbors_start = 108;
entryprobs = [0.99, 0.004, 0.001, 0.005] # [No entry, level1, level2, level3] - not taken from anything, just imposed.
entrants = [0, 1, 2, 3]
fields = 7; # if fields updated, update reshaping of state history
sim_start = 2;

#=
# TESTING --- To run a test, replace the value in the first bracket in mkt_fips = yearins[ ][1], and choose a year.

# Import Data
dataf = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
notmissing = findin(isna(dataf[:fipscode]), false);
dataf = dataf[notmissing, :];
yearins = [ [x; findfirst(dataf[:fipscode], x); findlast(dataf[:fipscode], x ); unique( dataf[findfirst(dataf[:fipscode], x):findlast(dataf[:fipscode], x ) , :year]  ) ] for x in unique(dataf[:fipscode])  ]
mkt_fips = yearins[145][1]
year = 2011
fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid]))

# There is a trick here in making sure this has the right size!  Need to compute the size of fids first.
state_history = [zeros(1, fields*size(fids)[1]) 1 0 0 0; zeros(T, fields*(size(fids)[1]) + 4)]

Simulator(dataf, year, mkt_fips, state_history, T = 100, sim_start = 2)
=#

function Simulator(dataf::DataFrame, year::Int64, mkt_fips::Int64,  state_history::Array{Float64,2}; T = 100, sim_start = 2)
  level1 = dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:level1_hospitals0][1]
  level2 = dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:level2solo_hospitals0][1]
  level3 = dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:level3_hospitals0][1]
  fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid]))
  for n in 1:size(fids)[1]
    el = fids[n]
    a = ((dataf[:,:fid].==el)&(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year))
    state_history[1, (n-1)*fields + 1] = el # Change
    state_history[1, (n-1)*fields + 2] = dataf[a,:act_int][1]
    state_history[1, (n-1)*fields + 3] = dataf[a,:act_solo][1]
    state_history[1, (n-1)*fields + 4] = 1 #probability is 1 for the first action
    state_history[1, (n-1)*fields + 5] = 10 # No action at the first period?  Or should it be 10?
    #state_history[1, (n-1)*fields + 6] = # whatever demand is
    state_history[1, (n-1)*fields + 7] = 0 #record the perturbation 0/1
  end
  # Record aggregte initial values
  state_history[1, (size(fids)[1])*fields+1] = level1 ;
  state_history[1, (size(fids)[1])*fields+2] = level2 ;
  state_history[1, (size(fids)[1])*fields+3] = level3 ;
  state_history[1, (size(fids)[1])*fields+4] = 1; # initial probability.

  # Start simulation here:
  for i = sim_start:T+1
    if i%50 == 0
        println(i)
    end
    fids = sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid])) # needs to be updated each round to catch entrants
        for fid in 1:size(fids)[1] # this has to be handled separately for each hospital, due to the geography issue
          el = fids[fid] # the dataframe is mutable.
          a = ((dataf[:,:fid].==el)&(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year))
          if sum(a) > 1
            println("two entries for ", el, " ", year, " ", mkt_fips)
          end
          if  ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (0,0)) # level 1, actions:
            probs1 = logitest((0,0), level1, level2, level3, convert(Array, [dataf[a,:lev105][1]; dataf[a,:lev205][1]; dataf[a,:lev305][1]; dataf[a,:lev1515][1]; dataf[a,:lev2515][1]; dataf[a,:lev3515][1]; dataf[a,:lev11525][1]; dataf[a,:lev21525][1]; dataf[a,:lev31525][1]]) )
            if probs1 == ValueException
              println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
              break
            end
            # all_hosp_probs[fid] = hcat(el, year_frame[a, :act_int], year_frame[a,:act_solo], 10, probs[1], 2, probs[2], 1, probs[3], 11, probs[4])
            # Reassign action choices -
            dataf[a, :choicenum0] = 10
            dataf[a, :pr_ch_0] = probs1[1]
            dataf[a, :choicenum1] = 2
            dataf[a, :pr_ch_1] = probs1[2]
            dataf[a, :choicenum2] = 1
            dataf[a, :pr_ch_2] = probs1[3]
            dataf[a, :choicenum3] = 11
            dataf[a, :pr_ch_3] = probs1[4]
            # Draw action:
            action1 = sample([10, 2, 1, 11] ,WeightVec([probs1[1], probs1[2], probs1[3], probs1[4]]))
            # Change things to reflect the action chosen:
            # record the probability of the action taken.
            if action1 == 10
              chprob = probs1[1]
              # no change to state in the aggregate
              # no change to distance-state of others
            elseif action1 == 2
              chprob = probs1[2]
              dataf[a,:act_int] = 1
              dataf[a,:act_solo] = 0
            elseif action1 == 1
              chprob = probs1[3]
              dataf[a,:act_int] = 0
              dataf[a,:act_solo] = 1
            elseif action1 == 11
              chprob = probs1[4]
              dataf[a,:act_int] = -999
              dataf[a,:act_solo] = -999
            else
              println("Fail")
              println("Bad Action Chosen by", el)
              break
            end
            # Set own distance counts to 0 for all categories
            (dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
            # This is the logical place to recompute the demand.
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields + 1] = el # Change
            state_history[i, (fid-1)*fields + 2] = dataf[a,:act_solo][1]
            state_history[i, (fid-1)*fields + 3] = dataf[a,:act_int][1]
            state_history[i, (fid-1)*fields + 4] = chprob
            state_history[i, (fid-1)*fields + 5] = action1
            #state_history[i, (fid-1)*fields + 6] = demand
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
          elseif ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (1,0)) #level 2, actions:
            # Can't evaluation as nexti here if they are initialized to negative 1
            probs2 = logitest((1,0), level1, level2, level3, convert(Array, [dataf[a,:lev105][1]; dataf[a,:lev205][1]; dataf[a,:lev305][1]; dataf[a,:lev1515][1]; dataf[a,:lev2515][1]; dataf[a,:lev3515][1]; dataf[a,:lev11525][1]; dataf[a,:lev21525][1]; dataf[a,:lev31525][1]]) )
            if probs2 == ValueException
              println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
              break
            end
            dataf[a, :choicenum0] = 5
            dataf[a, :pr_ch_0] = probs2[1]
            dataf[a, :choicenum1] = 10
            dataf[a, :pr_ch_1] = probs2[2]
            dataf[a, :choicenum2] = 6
            dataf[a, :pr_ch_2] = probs2[3]
            dataf[a, :choicenum3] = 11
            dataf[a, :pr_ch_3] = probs2[4]
            # Action:
            action2 = sample([5, 10, 6, 11] ,WeightVec([probs2[1], probs2[2], probs2[3], probs2[4]]))
            if action2 == 5
              chprob = probs2[1]
              dataf[a,:act_int] = 0
              dataf[a,:act_solo] = 0
            elseif action2 == 10
              chprob = probs2[2]
              # no change to state
            elseif action2 == 6
              chprob = probs2[3]
              dataf[a,:act_int] = 1
              dataf[a,:act_solo] = 0
            elseif action2 == 11
              chprob = probs2[4]
              dataf[a,:act_int] = -999
              dataf[a,:act_solo] = -999
            else
              println("Fail")
              println("Bad Action Chosen by", el)
              break
            end
            # Set own distance counts to 0 for all categories
            (dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields + 1] = el
            state_history[i, (fid-1)*fields + 2] = dataf[a,:act_solo][1]
            state_history[i, (fid-1)*fields + 3] = dataf[a,:act_int][1]
            state_history[i, (fid-1)*fields + 4] = chprob
            state_history[i, (fid-1)*fields + 5] = action2
            #state_history[i, (fid-1)*fields + 6] = demand
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
          elseif  ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (0,1)) #level 3, actions:
            probs3 = logitest((0,1), level1, level2, level3, convert(Array, [dataf[a,:lev105][1]; dataf[a,:lev205][1]; dataf[a,:lev305][1]; dataf[a,:lev1515][1]; dataf[a,:lev2515][1]; dataf[a,:lev3515][1]; dataf[a,:lev11525][1]; dataf[a,:lev21525][1]; dataf[a,:lev31525][1]]) )
            if probs3 == ValueException
              println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
              break
            end
            dataf[a, :choicenum0] = 4
            dataf[a, :pr_ch_0] = probs3[1]
            dataf[a, :choicenum1] = 3
            dataf[a, :pr_ch_1] = probs3[2]
            dataf[a, :choicenum2] = 10
            dataf[a, :pr_ch_2] = probs3[3]
            dataf[a, :choicenum3] = 11
            dataf[a, :pr_ch_3] = probs3[4]
            # Action:
            action3 = sample([4, 3, 10, 11] ,WeightVec([probs3[1], probs3[2], probs3[3], probs3[4]]))
            if action3 == 3
              chprob = probs3[2]
              dataf[a, :act_int] = 1
              dataf[a, :act_solo] = 0
            elseif action3 == 4
              chprob = probs3[1]
              dataf[a, :act_int] = 0
              dataf[a, :act_solo] =0
            elseif action3 == 10
              chprob = probs3[3]
              # no change to state
            elseif action3 == 11
              chprob = probs3[4]
              dataf[a, :act_int] = -999
              dataf[a, :act_solo] = -999
            else
              println("Fail")
              println("Bad Action Chosen by", el)
              break
            end
            # Set own distance counts to 0 for all categories
            (dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields + 1] = el
            state_history[i, (fid-1)*fields + 2] = dataf[a,:act_solo][1]
            state_history[i, (fid-1)*fields + 3] = dataf[a,:act_int][1]
            state_history[i, (fid-1)*fields + 4] = chprob
            state_history[i, (fid-1)*fields + 5] = action3
            #state_history[i, (fid-1)*fields + 6] = demand
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
          elseif ((dataf[a,:act_int][1], dataf[a,:act_solo][1]) == (-999,-999)) # has exited.
            # No new actions to compute, but record.
            state_history[i, (fid-1)*fields + 1] = el
            state_history[i, (fid-1)*fields + 2] = dataf[a,:act_solo][1]
            state_history[i, (fid-1)*fields + 3] = dataf[a,:act_int][1]
            state_history[i, (fid-1)*fields + 4] = 1 # exit is absorbing, so the choice prob is always 1
            state_history[i, (fid-1)*fields + 5] = 0 # no action is taken.
            #state_history[i, (fid-1)*fields + 6] = demand # no demand realized - exited.
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
            # Set own distance counts to 0 for all categories
            (dataf[a,:lev105], dataf[a,:lev205], dataf[a,:lev305], dataf[a,:lev1515], dataf[a,:lev2515], dataf[a,:lev3515], dataf[a,:lev11525], dataf[a,:lev21525], dataf[a,:lev31525]) = zeros(1,9)
          end
        end
        # Count facilities by distance and map results to neighboring hospitals
        # here the issue is that, for hospitals in neighboring counties, we haven't set the levX_YZ values to 0
        for f in fids
          own_fac = (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :act_solo][1], dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :act_int][1])
          for j = neighbors_start:(2):size(dataf)[2]
            if (!isna(dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j][1])) & (!isna(dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]))
              if (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), j+1][1] >0) & (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]< 5)
                if own_fac == (0,0)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev105] += 1
                elseif own_fac == (1,0)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev205] += 1
                elseif own_fac == (0,1)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev305] += 1
                elseif own_fac == (-999,-999)
                  # do nothing - firm exited.
                else
                  println("Bad Own Facility Code", f, j, own_fac)
                end
              elseif (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), j+1][1] >5) & (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]< 15)
                if own_fac == (0,0)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev1515] += 1
                elseif own_fac == (1,0)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev2515] += 1
                elseif own_fac == (0,1)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev3515] += 1
                elseif own_fac == (-999,-999)
                  # do nothing - firm exited.
                else
                  println("Bad Own Facility Code", f, j, own_fac)
                end
              elseif (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]> 15) & (dataf[(dataf[:fid].==f)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips),j+1][1]< 25)
                if own_fac == (0,0)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev11525] += 1
                elseif own_fac == (1,0)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev21525] += 1
                elseif own_fac == (0,1)
                  dataf[ (dataf[:fid].== j)&(dataf[:year].== year), :lev31525] += 1
                elseif own_fac == (-999,-999)
                  # do nothing - firm exited.
                else
                  println("Bad Own Facility Code", f, j, own_fac)
                end
              else
                println("Bad Distance at frame ", f, j)
              end
            end
          end
        end
    # Sum the levels for next period:
    level1 = 0; level2 = 0; level3 = 0;
    update_mkt = ((dataf[:year].==year)&(dataf[:fipscode].==mkt_fips)&(dataf[:act_int].!=-999)&(dataf[:act_solo].!=-999))
    total = sum(update_mkt)
    intens = sum(dataf[update_mkt, :act_int])
    solo = sum(dataf[update_mkt, :act_solo])
    nones = sum(update_mkt) - intens - solo
  #  println("Computing Market sizes")
        if (nones < 0) | (intens > total) | (solo > total) | (nones + intens + solo != total)
          println("Bad market size computations")
          println("Total ", total, " Level1 ", nones, " Level2 ", solo, " Level 3 ", intens, " Fips: ", mkt_fips, " Year ", year)
        else
            level1 = nones
            level2 = solo
            level3 = intens
        end

        # Entry draw
    entrypairs = hcat(entrants, entryprobs)
    entrantsp = WeightVec(entryprobs)
    newentrant = sample(entrants, entrantsp)
    entrantout = [newentrant, entrypairs[findfirst(entrypairs[:,1], newentrant), 2]]'
    b = ((dataf[:fipscode].== mkt_fips)&(dataf[:year].==year)); # Market-year observations, whether exited or not.
        if newentrant> 0
          println("Entry occurred")
          # Eventually fix the fact that the neighbors here are not going to be exactly right.
          # I need to add the fact that the entrants are not being recorded in the "neighbors" section so no one is counting distances to them.
          ent_lat = mean(dropna(dataf[b,:v15])) + rand(Normal(0, 0.1), 1) # 0.1 degrees latitude should be about 6-7 miles.
          ent_lon = mean(dropna(dataf[b,:v16])) + rand(Normal(0, 0.1), 1)
          newrow = dataf[b,:][1,:] # create new dataframe row, duplicating existing.  Takes first row of current
          newrow[:facility] = convert(UTF8String, "Entrant $mkt_fips $year")
          newrow[:fid] = sample(collect((maximum(dataf[:fid])+1):(maximum(dataf[:fid])+5))) # new facility will always have largest fid
          newrow[:id] = sample(collect((maximum(dataf[:id])+1):(maximum(dataf[:id])+5)))
          newrow[:fipscode] = mkt_fips
          newrow[:location] = convert(UTF8String, "entrant - see v15 v16")
          newrow[:city] = convert(UTF8String, "Entrant - unspecified")
          newrow[:firstyear] = year
          newrow[:v15] = ent_lat
          newrow[:v16] = ent_lon
          if newentrant == 1
            level1 += 1
            newrow[:act_int] = 0
            newrow[:act_solo] = 0
          elseif newentrant == 2
            level2 += 1
            newrow[:act_int] = 0
            newrow[:act_solo] = 1
          elseif newentrant == 3
            level3 += 1
            newrow[:act_int] = 0
            newrow[:act_solo] = 1
          end
          # Handle appending these entrants to the (neighbor, distance) section
          # need to check all of the other fids in the market-year (in b)
          println("Distances to new entrants")
          for row in eachrow(dataf[b,:])
            if  (distance(ent_lat[1], ent_lon[1], row[:v15], row[:v16]) < 25)
              count = 0 # only want to make an entry once - this is a dumb way
              for c in neighbors_start:(2):size(dataf)[2]
                if (isna(row[c]))&(count == 0)
                  println("doing something")
                  # to make changes I need to search for the row in the original DF matching these characteristics.
                  dataf[(dataf[:fid].==row[:fid])&(dataf[:id].==row[:id])&(dataf[:fipscode].==row[:fipscode]), c ]= newrow[:fid]
                  dataf[(dataf[:fid].==row[:fid])&(dataf[:id].==row[:id])&(dataf[:fipscode].==row[:fipscode]), c+1]= distance(ent_lat[1], ent_lon[1], row[:v15], row[:v16])
                  count += 1
                end
              end
            end
          end
          # Append the neighbors to the new entrant's frame too
          println("appending neighbors to new entrant's frame row")
          for elem in 1:size(fids)[1]
            neighb = fids[elem]
            neighb_lat = dataf[(dataf[:fid].==neighb)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :v15][1]
            neighb_lon = dataf[(dataf[:fid].==neighb)&(dataf[:year].==year)&(dataf[:fipscode].==mkt_fips), :v16][1]
            td = distance(ent_lat[1], ent_lon[1], neighb_lat, neighb_lon)
            if  td < 25
              newrow[neighbors_start+2*(elem-1)] = neighb # fid
              newrow[neighbors_start+2*(elem-1)+1] = td
              println("Success")
              if (td > 0) & (td < 5)
                if (newrow[:act_solo][1], newrow[:act_int][1]) == (0,0)
                  newrow[:lev105] += 1
                elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (1,0)
                  newrow[:lev205] += 1
                elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (0,1)
                  newrow[:lev305] += 1
                else
                  println("Bad facility in Entrant 1")
                end
              elseif (td > 5) & (td < 15)
                if (newrow[:act_solo][1], newrow[:act_int][1]) == (0,0)
                  newrow[:lev1515] += 1
                elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (1,0)
                  newrow[:lev2515] += 1
                elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (0,1)
                  newrow[:lev3515] += 1
                else
                  println("Bad facility in Entrant 2")
                end
              elseif (td > 15) & (td < 25)
                if (newrow[:act_solo][1], newrow[:act_int][1]) == (0,0)
                  newrow[:lev11525] += 1
                elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (1,0)
                  newrow[:lev21525] += 1
                elseif (newrow[:act_solo][1], newrow[:act_int][1]) == (0,1)
                  newrow[:lev31525] += 1
                else
                  println("Bad facility in Entrant 3")
                end
              else
                println("Bad distance measured from entrant")
              end
            end
          end
          println("appending entry")
          # Add the new record to the dataframe.
          append!(dataf, newrow)
          # append value to fids
          push!(fids, newrow[:fid][1])
          # Reshape state history: fid, solo state, int state, probability of choice,  action chosen, XXXX demand, perturbed. [newrow[:fid], 999, 999, 0, 1, 0]
          println("reshaping state history")
          # The problem is the size computation right here - figure it out.
          state_history = vcat(hcat(state_history[1:i,1:end-4], repmat([newrow[:fid][1] 999 999 1 0 0 0], i, 1), state_history[1:i, end-3:end]), zeros((T-i+1), size(fids)[1]*fields+4 ))
        end
    # Aggregate Probability of Action:
    tprob = 1
    for els in 4:fields:(size(state_history[i,:])[2]-4) # 4 is the relevant column
      if (state_history[i,els] <= 0) | (state_history[i,els]>1)
          println("Bad probability at row ", i, " ", els, " ", state_history[i,els-3 ], " prob ", state_history[i,els] )
          break
      else
        tprob = tprob*state_history[i,els]
      end
    end
    # Accounting for entry in aggregate prob of action -
    if newentrant >0
      if newentrant == 1
        tprob = tprob*entryprobs[2]
      elseif newentrant == 2
        tprob = tprob*entryprobs[3]
      elseif newentrant == 3
        tprob = tprob*entryprobs[4]
      end
    end
  #  println("computing aggregate transition probability")
    state_history[i, (size(fids)[1])*fields+4] = state_history[i-1, (size(fids)[1])*fields+4]*tprob  #prob of ending up at previous state * current transition prob

    # Total number of firms
    state_history[i, (size(fids)[1])*fields+1] = level1 ;
    state_history[i, (size(fids)[1])*fields+2] = level2 ;
    state_history[i, (size(fids)[1])*fields+3] = level3 ;

  end
  return state_history
end