function PerturbSimulator(data::Matrix, privatepeoplesub::Matrix, privatedemandmodelparameters::Array{Float64, 2}, medicaidpeoplesub::Matrix, medicaiddemandmodelparameters::Array{Float64, 2}, year::Int64, mkt_fips::Int64, pfid::Int64; entryprobs = [0.9895, 0.008, 0.0005, 0.002], entrants = [0, 1, 2, 3], disturb = 0.05,  T = 100, sim_start = 2, fields = 8, neighbors_start = 108, fipscodeloc = 78, yearloc = 75, level1_hospitals0loc = 11, fidloc = 74, level2solo_hospitals0loc = 10, level3_hospitals0loc = 9, act_intloc = 79, act_sololoc = 80, lev105loc = 97, lev205loc = 98, lev305loc = 99, lev1515loc = 101, lev2515loc = 102, lev3515loc = 103, lev11525loc = 105, lev21525loc = 106, lev31525loc = 107, v15loc = 94, v16loc = 95, idloc = 1, facilityloc = 82, locationloc = 88, firstyearloc = 91, cityloc = 85,  TotalBeds2loc = 20, TotalBeds1loc = 4)
  marketyear = (data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year) # It is a major speed up to do this once.
  level1 = data[marketyear,level1_hospitals0loc][1]
  level2 = data[marketyear,level2solo_hospitals0loc][1]
  level3 = data[marketyear,level3_hospitals0loc][1]
  fids = convert(Vector{Int64}, sort!(unique(data[marketyear,fidloc])))
  state_history = zeros(T+1, fields*size(fids,1)+4) #  0.000060 seconds (9 allocations: 58.688 KB)
  state_history[1, end-3] = 1 # set probability of initial outcome at 1

  for n in 1:size(fids,1)
    el = fids[n]
    a = (data[:,fidloc].==el)&marketyear
    state_history[1, (n-1)*fields+1] = el # Change
    state_history[1, (n-1)*fields+2] = data[a,act_intloc][1]
    state_history[1, (n-1)*fields+3] = data[a,act_sololoc][1]
    state_history[1, (n-1)*fields+4] = 1 #probability is 1 for the first action
    state_history[1, (n-1)*fields+5] = 10 # No action at the first period?  Or should it be 10?
    #state_history[1, (n-1)*fields + 6] = private # whatever demand is
    #state_history[1, (n-1)*fields + 7] = MEDICAID  # whatever demand is
    if el == pfid
      state_history[1, (n-1)*fields+8] = 1 #record the perturbation 0/1
    else
      state_history[1, (n-1)*fields+8] = 0 # not perturbed
    end
  end
  # Record aggregte initial values
  state_history[1, size(fids,1)*fields+1] = level1 ;
  state_history[1, size(fids,1)*fields+2] = level2 ;
  state_history[1, size(fids,1)*fields+3] = level3 ;
  state_history[1, size(fids,1)*fields+4] = 1; # initial probability.
  privatepeoplesub =rowchange(state_history[1,:], fids, privatepeoplesub) # fast - see DemandModel.jl. 10 fids: 0.148658 seconds (2.17 M allocations: 88.855 MB).  Probably improvable.
  medicaidpeoplesub =rowchange(state_history[1,:], fids, medicaidpeoplesub)
  emp_arr = Array{Float64, 2}()
  privaterealized_d = countmap(DemandModel(privatepeoplesub, privatedemandmodelparameters, emp_arr)) # maps chosen hospitals to counts.  Speed not great: 0.632453 seconds (757.94 k allocations: 347.089 MB, 36.85% gc time) / 10 hospitals.
  medicaidrealized_d = countmap(DemandModel(medicaidpeoplesub, medicaiddemandmodelparameters, emp_arr))
  for fid_i in 1:fields:size(state_history[1,:])[2]-4
    fid = state_history[1,fid_i]
    privatedemand_re =  try
      privaterealized_d[fid]
    catch y
      if isa(y, KeyError)
        privatedemand_re = 0 # write the demand out as -1 to keep track of failure to find val.
      else
        privatedemand_re = privaterealized_d[fid]
      end
    end
    state_history[1,fid_i+5] = privatedemand_re  #CHECK INDEX
    medicaiddemand_re =  try
      medicaidrealized_d[fid]
    catch y
      if isa(y, KeyError)
        medicaiddemand_re = 0 # write the demand out as -1 to keep track of failure to find val.
      else
        medicaiddemand_re = medicaidrealized_d[fid]
      end
    end
    state_history[1,fid_i+6] = medicaiddemand_re #CHECK INDEX
  end # Whole block: 0.003834 seconds (3.86 k allocations: 148.006 KB)
  # Start simulation here:
  total_entrants = Array{Float64, 2}() #must be defined outside the loop
  for i = sim_start:T+1
    marketyear = (data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year)
    fids = convert(Vector{Int64}, sort!(unique(data[marketyear,fidloc]))) # needs to be updated each round to catch entrants
        for fid in 1:size(fids,1) # this has to be handled separately for each hospital, due to the geography issue
          el = fids[fid] # the dataframe is mutable.
          a = ((data[:,fidloc].==el)&marketyear)
          prev_state = (state_history[i-1, (fid-1)*fields+2], state_history[i-1, (fid-1)*fields+3])
          if  prev_state ==  (0,0) # level 1, actions:
            probs1 = logitest((0,0), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
            if el == pfid
              probs1 = perturb(probs1, disturb, false) #perturb(probs::Array, eps::Float64, control::Bool)
            end
            # Draw action:
            action1 = sample([10, 2, 1, 11] ,WeightVec([probs1[1], probs1[2], probs1[3], probs1[4]]))
            # Change things to reflect the action chosen:
            # record the probability of the action taken.
            nextint = 0; nextsolo = 0; #reflects the current state
            if action1 == 10
              chprob = probs1[1]
              # no change to state in the aggregate
              # no change to distance-state of others
            elseif action1 == 2
              chprob = probs1[2]
              nextsolo = 0
              nextint = 1
            elseif action1 == 1
              chprob = probs1[3]
              nextsolo = 1
              nextint = 0
            elseif action1 == 11
              chprob = probs1[4]
              nextsolo = -999
              nextint = -999
            else
              # println("Fail")
              # println("Bad Action Chosen by", el)
              break
            end
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9)
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields+1] = el # Change
            state_history[i, (fid-1)*fields+2] = nextsolo
            state_history[i, (fid-1)*fields+3] = nextint
            state_history[i, (fid-1)*fields+4] = chprob
            state_history[i, (fid-1)*fields+5] = action1
            #state_history[i, (fid-1)*fields + 6] = demand
            #state_history[i, (fid-1)*fields + 7] = demand
            if el == pfid
              state_history[1, (fid-1)*fields+8] = 1 #record the perturbation 0/1
            else
              state_history[1, (fid-1)*fields+8] = 0 # not perturbed
            end
          elseif prev_state == (1,0) #level 2, actions:
            probs2 = logitest((1,0), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
            if el == pfid
              probs2 = perturb(probs2, disturb, false) #perturb(probs::Array, eps::Float64, control::Bool)
            end
            # Action:
            action2 = sample([5, 10, 6, 11] ,WeightVec([probs2[1], probs2[2], probs2[3], probs2[4]]))
            # Current state:
            nextint = 0; nextsolo = 1;
            if action2 == 5
              chprob = probs2[1]
              nextint = 0
              nextsolo = 0
            elseif action2 == 10
              chprob = probs2[2]
              # no change to state
            elseif action2 == 6
              chprob = probs2[3]
              nextint = 1
              nextsolo = 0
            elseif action2 == 11
              chprob = probs2[4]
              nextint = -999
              nextsolo = -999
            else
              # println("Fail")
              # println("Bad Action Chosen by", el)
              break
            end
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9)
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields+1] = el
            state_history[i, (fid-1)*fields+2] = nextsolo
            state_history[i, (fid-1)*fields+3] = nextint
            state_history[i, (fid-1)*fields+4] = chprob
            state_history[i, (fid-1)*fields+5] = action2
            #state_history[i, (fid-1)*fields + 6] = demand
            #state_history[i, (fid-1)*fields + 7] = demand
            if el == pfid
              state_history[1, (fid-1)*fields+8] = 1 #record the perturbation 0/1
            else
              state_history[1, (fid-1)*fields+8] = 0 # not perturbed
            end
          elseif  prev_state == (0,1) #level 3, actions:
            probs3 = logitest((0,1), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
            if el == pfid
              probs3 = perturb(probs3, disturb, false) #perturb(probs::Array, eps::Float64, control::Bool)
            end
            # Action:
            action3 = sample([4, 3, 10, 11] ,WeightVec([probs3[1], probs3[2], probs3[3], probs3[4]]))
            nextint = 1; nextsolo = 0; #current state
            if action3 == 3
              chprob = probs3[2]
              nextsolo = 1
              nextint = 0
            elseif action3 == 4
              chprob = probs3[1]
              nextint = 0
              nextsolo = 0
            elseif action3 == 10
              chprob = probs3[3]
              # no change to state
            elseif action3 == 11
              chprob = probs3[4]
              nextint = -999
              nextsolo = -999
            else
              # println("Fail")
              # println("Bad Action Chosen by", el)
              break
            end
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9)
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields+1] = el
            state_history[i, (fid-1)*fields+2] = nextsolo
            state_history[i, (fid-1)*fields+3] = nextint
            state_history[i, (fid-1)*fields+4] = chprob
            state_history[i, (fid-1)*fields+5] = action3
            #state_history[i, (fid-1)*fields + 6] = demand
            #state_history[i, (fid-1)*fields + 7] = demand
            if el == pfid
              state_history[1, (fid-1)*fields+8] = 1 #record the perturbation 0/1
            else
              state_history[1, (fid-1)*fields+8] = 0 # not perturbed
            end
          elseif prev_state == (-999,-999) # has exited.
            # No new actions to compute, but record.
            state_history[i, (fid-1)*fields+1] = el
            state_history[i, (fid-1)*fields+2] = -999
            state_history[i, (fid-1)*fields+3] = -999
            state_history[i, (fid-1)*fields+4] = 1 # exit is absorbing, so the choice prob is always 1
            state_history[i, (fid-1)*fields+5] = 0 # no action is taken.
            #state_history[i, (fid-1)*fields + 6] = demand # no demand realized - exited.
            #state_history[i, (fid-1)*fields + 7] = demand
            if el == pfid
              state_history[1, (fid-1)*fields+8] = 1 #record the perturbation 0/1
            else
              state_history[1, (fid-1)*fields+8] = 0 # not perturbed
            end
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9)
          end
        end
        # Count facilities by distance and map results to neighboring hospitals
        # here the issue is that, for hospitals in neighboring counties, we haven't set the levX_YZ values to 0
        #        print("part 2", "\n")
        for f in fids
          hospmktyear = (data[ :,fidloc].==f)&marketyear
          own_fac = (data[hospmktyear, act_sololoc][1], data[hospmktyear, act_intloc][1])
          for j = neighbors_start:(2):size(data, 2) # each row has appended to it a list of the hospital's neighbors.
            # starting at neighbors_start (col 108), we have (fid, distance) pairs in columns (or missing vals) all the way to the end of the dataframe
            nbhmktyear = (data[:,fidloc].==data[hospmktyear, j])&marketyear
            if sum(nbhmktyear)>0 # If the neighbor/market/year set isn't empty
              if (data[ nbhmktyear, act_intloc][1] != -999)&(data[ nbhmktyear, act_sololoc][1]!=-999) # if the neighbor hasn't exited
                if (!(data[hospmktyear, j][1] == 0)) & (!(data[hospmktyear,j+1][1] == 0)) # If the neighbor isn't missing, j - neighbor fid, j+1 - neighbor distance: if non-zero, proceed
                  if (data[hospmktyear, j+1][1] >0) & (data[hospmktyear,j+1][1]< 5) # distance to neighbor
                    if own_fac == (0,0)
                      data[ hospmktyear, lev105loc] += 1
                    elseif own_fac == (1,0)
                      data[ hospmktyear, lev205loc] += 1
                    elseif own_fac == (0,1)
                      data[ hospmktyear, lev305loc] += 1
                    elseif own_fac == (-999,-999)
                      # do nothing - firm exited.
                    else
                     #println("Bad Own Facility Code", f, j, own_fac)
                    end
                  elseif (data[hospmktyear, j+1][1] >5) & (data[hospmktyear,j+1][1]< 15)
                    if own_fac == (0,0)
                      data[ hospmktyear, lev1515loc] += 1
                    elseif own_fac == (1,0)
                      data[ hospmktyear, lev2515loc] += 1
                    elseif own_fac == (0,1)
                      data[ hospmktyear, lev3515loc] += 1
                    elseif own_fac == (-999,-999)
                      # do nothing - firm exited.
                    else
                      # println("Bad Own Facility Code", f, j, own_fac)
                    end
                  elseif (data[hospmktyear,j+1][1]> 15) & (data[hospmktyear,j+1][1]< 25)
                    if own_fac == (0,0)
                      data[ hospmktyear, lev11525loc] += 1
                    elseif own_fac == (1,0)
                      data[ hospmktyear, lev21525loc] += 1
                    elseif own_fac == (0,1)
                      data[ hospmktyear, lev31525loc] += 1
                    elseif own_fac == (-999,-999)
                      # do nothing - firm exited.
                    else
                      # println("Bad Own Facility Code", f, j, own_fac)
                    end
                  else
                    #  println("Bad Distance at frame ", f, j)
                  end
                end
              end
            end
          end
        end
    #    print("part 3", "\n")
            # Entry draw
    entrypairs = hcat(entrants, entryprobs)
    newentrant = sample(entrants, WeightVec(entryprobs))
    entrantout = [newentrant, entrypairs[findfirst(entrypairs[:,1], newentrant), 2]]'
    b = ((data[ :,fipscodeloc].== mkt_fips)&(data[ :,yearloc].==year)); # Market-year observations, whether exited or not.
        if newentrant> 0
          # I need to add the fact that the entrants are not being recorded in the "neighbors" section so no one is counting distances to them.
          ent_lat = mean((data[b,v15loc])) + rand(Normal(0, 0.1), 1) # 0.1 degrees latitude should be about 6-7 miles.
          ent_lon = mean((data[b,v16loc])) + rand(Normal(0, 0.1), 1)
          newrow = data[b,:][1,:] # create new row, duplicating existing.  Takes first row of current
          newrow[facilityloc] = convert(UTF8String, "Entrant $mkt_fips $year")
          newrow[fidloc] = maximum(data[:,fidloc])+5 # new facility will always have largest fid
          newrow[idloc] = -(maximum(data[:,idloc])+5) #
          newrow[fipscodeloc] = mkt_fips
          newrow[locationloc] = convert(UTF8String, "entrant - see v15 v16")
          newrow[cityloc] = convert(UTF8String, "Entrant - unspecified")
          newrow[firstyearloc] = year
          newrow[v15loc] = ent_lat
          newrow[v16loc] = ent_lon
          # Take the size as the mean bed number from neighboring hospitals.  There is no field for this in dataf, unfortunately.
# THIS IS WRONG
          entrantbeds = convert(Int, floor(mean( unique(unique(privatepeoplesub[ privatepeoplesub[:,TotalBeds1loc].>0 ,TotalBeds1loc])) ))) # 0.000021 seconds (28 allocations: 1.406 KB)
          if newentrant == 1
            level1 += 1
            newrow[act_intloc] = 0
            newrow[act_sololoc] = 0
          elseif newentrant == 2
            level2 += 1
            newrow[act_intloc] = 0
            newrow[act_sololoc] = 1
          elseif newentrant == 3
            level3 += 1
            newrow[act_intloc] = 0
            newrow[act_sololoc] = 1
          end
          # Define the data needed for an entrant to compute the demand model for the individuals.  6 components.
          #        print("part 3.25 \n")
          entrant_data = [newrow[fidloc] newrow[act_intloc] newrow[act_sololoc] entrantbeds ent_lat ent_lon]
          if maximum(size(total_entrants))<=1
            total_entrants = entrant_data
          else
            total_entrants = [total_entrants entrant_data]
          end
          # Handle appending these entrants to the (neighbor, distance) section
          # need to check all of the other fids in the market-year (in b)
          for row in findfirst(marketyear, 1):findlast(marketyear, 1) # find the first and last row in the market year - this will return the row number in the whole array
            if  (distance(ent_lat[1], ent_lon[1], data[row,v15loc][1], data[row,v16loc][1]) < 25)
              valfn = findfirst(data[row, neighbors_start-1:2:end], 0) # take the market via data[b,:], go through it via "row", but the number returned needs to be multiplied by 2 and added to neighbors_start!
              if (valfn > 0)&(valfn < 36)
                nbloc = neighbors_start + 2*valfn
                if nbloc < size(data,2)
                  data[row, nbloc]= newrow[fidloc] # potentially an unknown error occurs here sometimes.
                  data[row, nbloc+1]= distance(ent_lat[1], ent_lon[1], newrow[v15loc][1], newrow[v16loc][1])
                else
                  println("There's an error in the size of nbloc, line 300 ", nbloc)
                  println("size data ", size(data))
                  println("period? ", i)
                end
               end
            end
          end
          # Append the neighbors to the new entrant's frame too
          #  println("appending neighbors to new entrant's frame row")
          for elem in 1:size(fids,1)
            neighb = fids[elem]
            neighb_row = (data[ :,fidloc].==neighb)&marketyear
            neighb_lat = data[neighb_row, v15loc][1]
            neighb_lon = data[neighb_row, v16loc][1]
            td = distance(ent_lat[1], ent_lon[1], neighb_lat[1], neighb_lon[1])
            #    print("computed distance: ", td, "\n")
            if  td < 25
              # A strange problem could arise here if in some simulation there are more than 43 or so entrants - there won't be enough room in the row.  Not a huge worry.
              newrow[neighbors_start+2*(elem-1)] = neighb # fid
              newrow[neighbors_start+2*(elem-1)+1] = td
              entrant_state = (newrow[act_sololoc][1], newrow[act_intloc][1])
              if (td > 0) & (td < 5)
                if entrant_state == (0,0)
                  newrow[lev105loc] += 1
                elseif entrant_state == (1,0)
                  newrow[lev205loc] += 1
                elseif entrant_state == (0,1)
                  newrow[lev305loc] += 1
                else
                  #              println("Bad facility in Entrant 1")
                end
              elseif (td > 5) & (td < 15)
                if entrant_state == (0,0)
                  newrow[lev1515loc] += 1
                elseif entrant_state == (1,0)
                  newrow[lev2515loc] += 1
                elseif entrant_state == (0,1)
                  newrow[lev3515loc] += 1
                else
                  #            println("Bad facility in Entrant 2")
                end
              elseif (td > 15) & (td < 25)
                if entrant_state == (0,0)
                  newrow[lev11525loc] += 1
                elseif entrant_state == (1,0)
                  newrow[lev21525loc] += 1
                elseif entrant_state == (0,1)
                  newrow[lev31525loc] += 1
                else
                  #            println("Bad facility in Entrant 3")
                end
              else
                #          println("Bad distance measured from entrant")
              end
            end
          end
          #    println("appending entry")
          # Add the new record to the dataframe.
          #      print("part 3.5 \n")
          data = [data; newrow]
          # append value to fids
          push!(fids, newrow[fidloc][1])
          # Reshape state history: fid, solo state, int state, probability of choice,  action chosen, XXXX demand, perturbed. [newrow[:fid], 999, 999, 0, 1, 0]
          #      println("reshaping state history")
          # The problem is the size computation right here - figure it out.
          state_history = vcat( hcat(state_history[1:i-1, 1:end-4], repmat([newrow[fidloc][1] 999 999 1 0 0 0 0], i-1, 1), state_history[1:i-1,end-3:end]), hcat(state_history[i,1:end-4], [newrow[fidloc] newrow[act_sololoc] newrow[act_intloc] entrantout[2] 0 0 0 0], state_history[i, end-3:end]) ,zeros((T-i+1), size(fids,1)*fields+4) ) # 0.000126 seconds (80 allocations: 132.234 KB) ??
        end
    #    print("part 4", "\n")
    # Aggregate Probability of Action:
    tprob = 1
    #  print("Part 4", "\n")
      for els in 4:fields:(size(state_history[i,:], 2)-4) # 4 is the relevant column
        if (state_history[i,els] <= 0) | (state_history[i,els]>1)
            break
        else
          tprob = tprob*state_history[i,els]
        end
      end
    # Accounting for entry in aggregate prob of action -
        if newentrant == 0
          tprob = tprob*entryprobs[1]
        elseif newentrant == 1
          tprob = tprob*entryprobs[2]
        elseif newentrant == 2
          tprob = tprob*entryprobs[3]
        elseif newentrant == 3
          tprob = tprob*entryprobs[4]
        end
        # Sum the levels for next period:
      #  print("part 5", "\n")
      level1 = 0; level2 = 0; level3 = 0;
      update_mkt = state_history[i, 1:fields:end-4]
      total = maximum(size(update_mkt))
      solo = sum(state_history[i, 2:fields:end-4].==1) # counts solo facilities from records in state history
      intens = sum(state_history[i, 3:fields:end-4].==1) # counts intensives from state history directly
      nones = sum((state_history[i, 2:fields:end-4].==0)&(state_history[i, 3:fields:end-4].==0))
      level1 = nones
      level2 = solo
      level3 = intens
    # Here is the place to do demand - map the results above out to the demand model
    #=
    - Write the fids out to an array
    - Statehistory row has been computed now - can use that. Current market state given by statehistory[i,:]
    - Rows should be changed by rowchange(staterow::Array{Float64,2}, choicerow::DataFrame) for every row in peoplesub
    - At the end call DemandModel(people::DataFrame, modelparameters::Array{Float64, 2}) on the result.
    - Obtain demand and map it into state_history
    =#
    privatepeoplesub =rowchange(state_history[i-1,:], fids, privatepeoplesub) # fast - see DemandModel.jl. 10 fids: 0.148658 seconds (2.17 M allocations: 88.855 MB).  Probably improvable.
    medicaidpeoplesub =rowchange(state_history[i-1,:], fids, medicaidpeoplesub)
    privaterealized_d = countmap(DemandModel(privatepeoplesub, privatedemandmodelparameters, total_entrants)) # maps chosen hospitals to counts.  Speed not great: 0.632453 seconds (757.94 k allocations: 347.089 MB, 36.85% gc time) / 10 hospitals.
    medicaidrealized_d = countmap(DemandModel(medicaidpeoplesub, medicaiddemandmodelparameters, total_entrants))
    for fid_i in 1:fields:size(state_history[i,:])[2]-4
      fid = state_history[i,fid_i]
      privatedemand_re =  try
        privaterealized_d[fid]
      catch y
        if isa(y, KeyError)
          privatedemand_re = 0 # write the demand out as -1 to keep track of failure to find val.
        else
          privatedemand_re = privaterealized_d[fid]
        end
      end
      state_history[i,fid_i+5] = privatedemand_re  #CHECK INDEX
      medicaiddemand_re =  try
        medicaidrealized_d[fid]
      catch y
        if isa(y, KeyError)
          medicaiddemand_re = 0 # write the demand out as -1 to keep track of failure to find val.
        else
          medicaiddemand_re = medicaidrealized_d[fid]
        end
      end
      state_history[i,fid_i+6] = medicaiddemand_re #CHECK INDEX
    end
    #  println("computing aggregate transition probability")
    state_history[i, (size(fids)[1])*fields+4] = state_history[i-1, (size(fids)[1])*fields+4]*tprob  #prob of ending up at previous state * current transition prob
    # Total number of firms
    state_history[i, (size(fids)[1])*fields+1] = level1 ;
    state_history[i, (size(fids)[1])*fields+2] = level2 ;
    state_history[i, (size(fids)[1])*fields+3] = level3 ;
  end
  privatepeoplesub =rowchange(state_history[1,:], fids, privatepeoplesub) # fast - see DemandModel.jl. 10 fids: 0.148658 seconds (2.17 M allocations: 88.855 MB).  Probably improvable.
  medicaidpeoplesub =rowchange(state_history[1,:], fids, medicaidpeoplesub) # reset the people to the original state
  index = findfirst(state_history[1,:], pfid) # only return the state history of the perturbed guy.
  subset = hcat(state_history[:, index:index+7], state_history[:, end-3:end])
  return subset # state_history # change order of subset and state_history to return the whole perturbed state.
end






#=
# TESTING --- To run a test, replace the value in the first bracket in mkt_fips = yearins[ ][1], and choose a year.

# Import Data

yearins = [ [x; findfirst(data[:,fipscodeloc], x); findlast(data[:,fipscodeloc], x ); unique( data[findfirst(data[:,fipscodeloc], x):findlast(data[:,fipscodeloc], x ) , yearloc]  ) ] for x in unique(data[:, fipscodeloc])  ]
mkt_fips = yearins[125][1]
year = 2011
fids = sort!(unique(data[(data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year),fidloc]))
pfid = fids[1]
disturb = 0.05


state_history = PerturbSimulator(data, peoples, year, mkt_fips, demandmodelparameters, pfid, disturb = 0.05, T = 10, sim_start = 2)
# Find the right people:
peoplesub = people[fidfinder(convert(Array, fids)', people, "people"),:]

# To reset for repeated simulations:: (This eliminates entrants, all of which have negative id's)
dataf = dataf[dataf[:id].>= 0, :]
=#

#entryprobs = [0.9895, 0.008, 0.0005, 0.002]; entrants = [0, 1, 2, 3]; disturb = 0.05; T = 100; sim_start = 2; fields = 7; neighbors_start = 108; fipscodeloc = 78; yearloc = 75; level1_hospitals0loc = 11; fidloc = 74; level2solo_hospitals0loc = 10; level3_hospitals0loc = 9; act_intloc = 79; act_sololoc = 80; lev105loc = 97; lev205loc = 98; lev305loc = 99; lev1515loc = 101; lev2515loc = 102; lev3515loc = 103; lev11525loc = 105; lev21525loc = 106; lev31525loc = 107; v15loc = 94; v16loc = 95; idloc = 1; facilityloc = 82; locationloc = 88; firstyearloc = 91; cityloc = 85; TotalBeds2loc = 20; TotalBeds1loc = 4







###
