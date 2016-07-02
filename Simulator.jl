function Simulator(data::Matrix, peoplesub::Matrix, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2}; T = 100, sim_start = 2, fields = 7, neighbors_start = 108, entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002], fipscodeloc = 78, yearloc = 75, level1_hospitals0loc = 11, fidloc = 74, level2solo_hospitals0loc = 10, level3_hospitals0loc = 9, act_intloc = 79, act_sololoc = 80, lev105loc = 97, lev205loc = 98, lev305loc = 99, lev1515loc = 101, lev2515loc = 102, lev3515loc = 103, lev11525loc = 105, lev21525loc = 106, lev31525loc = 107, v15loc = 94, v16loc = 95, idloc = 1, facilityloc = 82, locationloc = 88, firstyearloc = 91, cityloc = 85,  TotalBeds2loc = 20, TotalBeds1loc = 4)
# think about changing the "data" matrix to be type Array{Float64, 2} - can drop strings in the Reboot.jl - see if speeds up.
  marketyear = (data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year) # It is a major speed up to do this once. (14.00 k allocations: 377.698 KB)
  level1 = data[marketyear,level1_hospitals0loc][1] # 7 allocations: 464 bytes
  level2 = data[marketyear,level2solo_hospitals0loc][1] # 8 allocations: 512 bytes
  level3 = data[marketyear,level3_hospitals0loc][1] # 8 allocations: 512 bytes
  fids = convert(Vector{Int64}, sort!(unique(data[marketyear,fidloc]))) # 0.006925 seconds 4.06 k allocations: 208.760 KB - slow.
  state_history = zeros(T+1, fields*size(fids,1) + 4) #  0.000060 seconds (9 allocations: 58.688 KB)
  state_history[1, end-3] = 1 # set probability of initial outcome at 1


  # Writes the values to the first row of the state history
  for n in 1:size(fids,1)
    el = fids[n]
    a = ((data[:,fidloc].==el)&marketyear) # 0.000446 seconds (6.93 k allocations: 168.484 KB)
    state_history[1, (n-1)*fields + 1] = el # Change
    state_history[1, (n-1)*fields + 2] = data[a,act_intloc][1] # 0.000020 seconds (7 allocations: 336 bytes)
    state_history[1, (n-1)*fields + 3] = data[a,act_sololoc][1] # 0.000024 seconds (7 allocations: 336 bytes)
    state_history[1, (n-1)*fields + 4] = 1 #probability is 1 for the first action
    state_history[1, (n-1)*fields + 5] = 10 # No action at the first period?  Or should it be 10?
    #state_history[1, (n-1)*fields + 6] = #  demand is not recorded here - recorded below.
    # state_history[1, (n-1)*fields + 7] = 0 #record the perturbation 0/1 - always 0 in equilibrium simulation.
  end
  # Record aggregte initial values
  state_history[1, size(fids,1)*fields+1] = level1 ; # 0.000006 seconds (4 allocations: 160 bytes)
  state_history[1, size(fids,1)*fields+2] = level2 ;
  state_history[1, size(fids,1)*fields+3] = level3 ;
  state_history[1, size(fids,1)*fields+4] = 1; # initial probability.
  # Compute initial demand here:
  peoplesub =rowchange(state_history[1,:], fids, peoplesub) # fast - see DemandModel.jl. 10 fids: 0.148658 seconds (2.17 M allocations: 88.855 MB).  Probably improvable.
  #DemandModel(people::DataFrame, frname::ASCIIString, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; maxfid = 11, ent_length = 6 )
  emp_arr = Array{Float64, 2}()
  realized_d = countmap(DemandModel(peoplesub, demandmodelparameters, emp_arr)) # maps chosen hospitals to counts.  Speed not great: 0.632453 seconds (757.94 k allocations: 347.089 MB, 36.85% gc time) / 10 hospitals.
  for fid_i in 1:fields:size(state_history[1,:])[2]-4
    fid = state_history[1,fid_i]
    demand_re =  try
      realized_d[fid]
    catch y
      if isa(y, KeyError)
        demand_re = 0 # write the demand out as -1 to keep track of failure to find val.
      else
        demand_re = realized_d[fid]
      end
    end
    state_history[1,fid_i+5] = demand_re
  end # Whole block: 0.003834 seconds (3.86 k allocations: 148.006 KB)
  # Start simulation here:
  total_entrants = Array{Float64, 2}() # define outside the loop, only for the first round.
  for i = sim_start:T+1
    marketyear = (data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year) # 0.000962 seconds (13.85 k allocations: 335.766 KB)
      fids = convert(Vector{Int64}, sort!(unique(data[marketyear,fidloc]))) # needs to be updated each round to catch entrants
        for fid in 1:size(fids,1) # this has to be handled separately for each hospital, due to the geography issue
          el = fids[fid] # takes the fid
          a = ((data[:,fidloc].==el)&marketyear)
          prev_state = (state_history[i-1, (fid-1)*fields + 2], state_history[i-1, (fid-1)*fields + 3]) #0.000012 seconds (8 allocations: 272 bytes)
          if prev_state ==  (0,0)                     # ((data[a,act_intloc][1], data[a,act_sololoc][1]) == (0,0)) # level 1, actions:
            probs1 = logitest((0,0), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] ) # 0.000200 seconds (419 allocations: 25.500 KB)
            if probs1 == ProjectModule.ValueException
      #        println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
              break
            end
            # Draw action:
            action1 = sample([10, 2, 1, 11] ,WeightVec([probs1[1], probs1[2], probs1[3], probs1[4]])) # @time 0.000020 seconds (11 allocations: 448 bytes)
            # Change things to reflect the action chosen:
            nextint = 0; nextsolo = 0; #reflects the current state
            if action1 == 10
              chprob = probs1[1]
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
  #            println("Fail")
  #            println("Bad Action Chosen by", el)
              break
            end # 0.000004 seconds (5 allocations: 176 bytes)
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9) #0.000076 seconds (53 allocations: 1.781 KB)
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields + 1] = el # writes fid
            state_history[i, (fid-1)*fields + 2] = nextsolo
            state_history[i, (fid-1)*fields + 3] = nextint
            state_history[i, (fid-1)*fields + 4] = chprob
            state_history[i, (fid-1)*fields + 5] = action1
            #state_history[i, (fid-1)*fields + 6] = demand # handled below, not here.
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
          elseif prev_state == (1,0)  #level 2, actions:
            # Can't evaluation as nexti here if they are initialized to negative 1
            probs2 = logitest((1,0), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
            if probs2 == ProjectModule.ValueException
  #            println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
              break
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
            state_history[i, (fid-1)*fields + 1] = el
            state_history[i, (fid-1)*fields + 2] = nextsolo
            state_history[i, (fid-1)*fields + 3] = nextint
            state_history[i, (fid-1)*fields + 4] = chprob
            state_history[i, (fid-1)*fields + 5] = action2
            #state_history[i, (fid-1)*fields + 6] = demand # Not recorded here - recorded below.
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
          elseif prev_state == (0,1)   #level 3, actions:
            probs3 = logitest((0,1), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
            if probs3 == ProjectModule.ValueException
    #          println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
              break
            end
            # Action:
            action3 = sample([4, 3, 10, 11] ,WeightVec([probs3[1], probs3[2], probs3[3], probs3[4]]))
            nextint = 1; nextsolo = 0; #current state
            if action3 == 3 # go to 2
              chprob = probs3[2]
              nextsolo = 1
              nextint = 0
            elseif action3 == 4 #go to 1
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
  #            println("Fail")
  #            println("Bad Action Chosen by", el)
              break
            end
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9)
            # write out state values - in blocks:
            state_history[i, (fid-1)*fields + 1] = el
            state_history[i, (fid-1)*fields + 2] = nextsolo
            state_history[i, (fid-1)*fields + 3] = nextint
            state_history[i, (fid-1)*fields + 4] = chprob
            state_history[i, (fid-1)*fields + 5] = action3
            #state_history[i, (fid-1)*fields + 6] = demand # Not recorded here - recorded below.
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
          elseif prev_state == (-999,-999)   # has exited.
            # No new actions to compute, but record.
  #          print("Exit recorded \n")
            state_history[i, (fid-1)*fields + 1] = el
            state_history[i, (fid-1)*fields + 2] = -999
            state_history[i, (fid-1)*fields + 3] = -999
            state_history[i, (fid-1)*fields + 4] = 1 # exit is absorbing, so the choice prob is always 1
            state_history[i, (fid-1)*fields + 5] = 0 # no action is taken.
            #state_history[i, (fid-1)*fields + 6] = demand # no demand realized - exited. Not recorded here - recorded below.
            state_history[i, (fid-1)*fields + 7] = 0 #perturbation
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9)
          end
        end # The whole block above took (for 10 facilities) 0.084701 seconds (125.78 k allocations: 4.188 MB)
        # Measure distances to neighbors
        #### Speed up opportunity.
        for f in fids # 0.159000 seconds (82.68 k allocations: 23.804 MB) - for 10 fids.  This can be sped up.
          hospmktyear = (data[ :,fidloc].==f)&marketyear # 0.000447 seconds (6.93 k allocations: 168.531 KB)
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
      #                println("Bad Own Facility Code", f, j, own_fac)
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
  #                    println("Bad Own Facility Code", f, j, own_fac)
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
    #                  println("Bad Own Facility Code", f, j, own_fac)
                    end
                  else
    #                println("Bad Distance at frame ", f, j)
                  end
                end
              end
            end
          end
        end


        # Entry draw
#        print("part 3 \n")
      entrypairs = hcat(entrants, entryprobs)
      newentrant = sample(entrants, WeightVec(entryprobs))
      entrantout = [newentrant, entrypairs[findfirst(entrypairs[:,1], newentrant), 2]]' # 0.000018 seconds (13 allocations: 560 bytes)
      b = ((data[ :,fipscodeloc].== mkt_fips)&(data[ :,yearloc].==year)); # Market-year observations, whether exited or not.
        if newentrant> 0
          # I need to add the fact that the entrants are not being recorded in the "neighbors" section so no one is counting distances to them.
          ent_lat = mean((data[b,v15loc])) + rand(Normal(0, 0.1), 1) # 0.000033 seconds (22 allocations: 848 bytes) 0.1 degrees latitude should be about 6-7 miles.
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
          # Just compute the mean over all beds in the state?  This needs to be fixed later.
          entrantbeds = convert(Int, floor(mean( unique(unique(peoplesub[ peoplesub[TotalBeds1loc].>0 ,TotalBeds1loc])) ))) # 0.000021 seconds (28 allocations: 1.406 KB)
          # This part IS necessary
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
          entrant_data = [newrow[fidloc] newrow[act_intloc] newrow[act_sololoc] entrantbeds ent_lat ent_lon]  # 0.000080 seconds (105 allocations: 4.531 KB)
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
          # println("appending neighbors to new entrant's frame row")
          for elem in 1:size(fids,1)
            neighb = fids[elem]
            neighb_row = (data[ :,fidloc].==neighb)&marketyear
            neighb_lat = data[neighb_row, v15loc][1]
            neighb_lon = data[neighb_row, v16loc][1]
            td = distance(ent_lat[1], ent_lon[1], neighb_lat[1], neighb_lon[1])
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
          # println("appending entry")
          # Add the new record to the dataframe.
          data = [data; newrow] # quick but costly: 0.007540 seconds (17 allocations: 9.541 MB)
          # append value to fids
          push!(fids, newrow[fidloc][1])
          # Reshape state history: fid, solo state, int state, probability of choice,  action chosen, XXXX demand, perturbed. [newrow[:fid], 999, 999, 0, 1, 0]
          #    println("reshaping state history")
          # The problem is the size computation right here - figure it out.  I'm not resizing this right.  The entry condition is not correct now.
          # OLD:  vcat(hcat(state_history[1:i,1:end-4], repmat([newrow[fidloc][1] 999 999 1 0 0 0], i, 1), state_history[1:i, end-3:end]), zeros((T-i+1), size(fids,1)*fields+4 ))
          # Want: vcat( hcat(state_history[1:i-1, 1:end-4], repmat([newrow[fidloc][1] 999 999 1 0 0 0], i-1, 1), state_history[1:i-1,end-3:end]), hcat(state_history[i,1:end-4], [newrow[fidloc] newrow[act_sololoc] newrow[act_intloc] entrantout[2] 0 0 0], state_history[i, end-3:end]) ,zeros((T-i+1), (size(fids,1)+1)*fields+4) ))
          state_history = vcat( hcat(state_history[1:i-1, 1:end-4], repmat([newrow[fidloc][1] 999 999 1 0 0 0], i-1, 1), state_history[1:i-1,end-3:end]), hcat(state_history[i,1:end-4], [newrow[fidloc] newrow[act_sololoc] newrow[act_intloc] entrantout[2] 0 0 0], state_history[i, end-3:end]) ,zeros((T-i+1), size(fids,1)*fields+4) ) # 0.000126 seconds (80 allocations: 132.234 KB) ??
# The above is super inefficient. Can we allocate a new matrix and set certain elements equal to the old ones?

        end
        # Aggregate Probability of Action:
      tprob = 1
#      print("Part 4", "\n")
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
      level1 = 0; level2 = 0; level3 = 0;
      update_mkt = ((data[:,yearloc].==year)&(data[:,fipscodeloc].==mkt_fips)&(data[:,act_intloc].!=-999)&(data[:,act_sololoc].!=-999)) #compute again for entrants/exiters
      total = sum(update_mkt)
      intens = sum(data[update_mkt, act_intloc])
      solo = sum(data[update_mkt, act_sololoc])
      nones = sum(update_mkt) - intens - solo
      #  println("Computing Market sizes")
        if (nones < 0) | (intens > total) | (solo > total) | (nones + intens + solo != total)
  #        println("Bad market size computations")
  #        println("Total ", total, " Level1 ", nones, " Level2 ", solo, " Level 3 ", intens, " Fips: ", mkt_fips, " Year ", year)
        else
            level1 = nones
            level2 = solo
            level3 = intens
        end
    #    print("Part 4.5 \n")
      # Here is the place to do demand - map the results above out to the demand model
      #=
      - Write the fids out to an array
      - Statehistory row has been computed now - can use that. Current market state given by statehistory[i,:]
      - Rows should be changed by rowchange(staterow::Array{Float64,2}, choicerow::DataFrame) for every row in peoplesub
      - At the end call DemandModel(people::DataFrame, modelparameters::Array{Float64, 2}) on the result.
      - Obtain demand and map it into state_history
      =#
      peoplesub = rowchange(state_history[i-1,:], fids, peoplesub)
      realized_d = countmap(DemandModel(peoplesub, demandmodelparameters, total_entrants)) # maps chosen hospitals to counts.
      for fid_i in 1:fields:size(state_history[i,:], 2)-4
        fid = state_history[i,fid_i]
        demand_re =  try
          realized_d[fid]
        catch y
          if isa(y, KeyError)
            demand_re = 0 # write the demand out as -1 to keep track of failure to find val.
          else
            demand_re = realized_d[fid]
          end
        end
        state_history[i,fid_i+5] = demand_re
      end
#      print("Part 5 \n")
      #  println("computing aggregate transition probability")
      state_history[i, (size(fids)[1])*fields+4] = state_history[i-1, (size(fids)[1])*fields+4]*tprob  #prob of ending up at previous state * current transition prob

      # Total number of firms
      state_history[i, (size(fids)[1])*fields+1] = level1 ;
      state_history[i, (size(fids)[1])*fields+2] = level2 ;
      state_history[i, (size(fids)[1])*fields+3] = level3 ;

    end
  peoplesub =rowchange(state_history[1,:], fids, peoplesub) # reset the people to the original state
  return state_history
end



#=
# TESTING --- To run a test, replace the value in the first bracket in mkt_fips = yearins[ ][1], and choose a year.


# Import Data

yearins = [ [x; findfirst(data[:,fipscodeloc], x); findlast(data[:,fipscodeloc], x ); unique( data[findfirst(data[:,fipscodeloc], x):findlast(data[:,fipscodeloc], x ) , yearloc]  ) ] for x in unique(data[:,fipscodeloc])  ]
mkt_fips = yearins[11][1]
year = 2005
fids = convert(Vector, sort!(unique(data[(data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year),fidloc])))



# To reset for repeated simulations:: (This eliminates entrants, all of which have negative id's)
data = data[data[:,idloc].>= 0, :]

Simulator(dataf, peoplesub, "peoplesub", year, mkt_fips, demandmodelparameters)


fipscodeloc = dataf.colindex.lookup[:fipscode]
fipscodeloc = 78;

yearloc = dataf.colindex.lookup[:year]
yearloc = 75;

level1_hospitals0loc = dataf.colindex.lookup[:level1_hospitals0]
level1_hospitals0loc = 11;
level2solo_hospitals0loc = dataf.colindex.lookup[:level2solo_hospitals0]
level2solo_hospitals0loc = 10
level3_hospitals0loc = dataf.colindex.lookup[:level3_hospitals0]
level3_hospitals0loc = 9
fidloc = dataf.colindex.lookup[:fid]
fidloc = 74
act_intloc = dataf.colindex.lookup[:act_int]
act_intloc = 79
act_sololoc = dataf.colindex.lookup[:act_solo]
act_sololoc = 80
lev105loc = dataf.colindex.lookup[:lev105]
lev105loc = 97
lev205loc = dataf.colindex.lookup[:lev205]
lev205loc = 98
lev305loc = dataf.colindex.lookup[:lev305]
lev305loc = 99
lev1515loc = dataf.colindex.lookup[:lev1515]
lev1515loc = 101
lev2515loc = dataf.colindex.lookup[:lev2515]
lev2515loc = 102
lev3515loc = dataf.colindex.lookup[:lev3515]
lev3515loc = 103
lev11525loc = dataf.colindex.lookup[:lev11525]
lev11525loc = 105
lev21525loc = dataf.colindex.lookup[:lev21525]
lev21525loc = 106
lev31525loc = dataf.colindex.lookup[:lev31525]
lev31525loc = 107
v15loc = dataf.colindex.lookup[:v15]
v15loc = 94
v16loc = dataf.colindex.lookup[:v16]
v16loc = 95
facilityloc = dataf.colindex.lookup[:facility]
facilityloc = 82
idloc = dataf.colindex.lookup[:id]
idloc = 1
locationloc = dataf.colindex.lookup[:location]
locationloc = 88
cityloc = dataf.colindex.lookup[:city]
cityloc = 85
firstyearloc = dataf.colindex.lookup[:firstyear]
firstyearloc = 91
TotalBeds1loc = people.colindex.lookup[:TotalBeds1]
TotalBeds1loc = 4
TotalBeds2loc = people.colindex.lookup[:TotalBeds2]
TotalBeds2loc = 20

T = 100; sim_start = 2; fields = 7; neighbors_start = 108; entrants = [0, 1, 2, 3]; entryprobs = [0.9895, 0.008, 0.0005, 0.002]
fipscodeloc = 78; yearloc = 75; level1_hospitals0loc = 11; TotalBeds2loc = 20; fidloc = 74; level2solo_hospitals0loc = 10; level3_hospitals0loc = 9; act_intloc = 79; act_sololoc = 80; lev105loc = 97; lev205loc = 98; lev305loc = 99; lev1515loc = 101; lev2515loc = 102; lev3515loc = 103; lev11525loc = 105; lev21525loc = 106; lev31525loc = 107; v15loc = 94; v16loc = 95; idloc = 1; facilityloc = 82; locationloc = 88; firstyearloc = 91; cityloc = 85;  TotalBeds2loc = 20; TotalBeds1loc = 4;
fipscodeloc = 78, yearloc = 75, level1_hospitals0loc = 11, TotalBeds2loc = 20, fidloc = 74, level2solo_hospitals0loc = 10, level3_hospitals0loc = 9, act_intloc = 79, act_sololoc = 80, lev105loc = 97, lev205loc = 98, lev305loc = 99, lev1515loc = 101, lev2515loc = 102, lev3515loc = 103, lev11525loc = 105, lev21525loc = 106, lev31525loc = 107, v15loc = 94, v16loc = 95, idloc = 1, facilityloc = 82, locationloc = 88, firstyearloc = 91, cityloc = 85,  TotalBeds2loc = 20, TotalBeds1loc = 4,
=#
