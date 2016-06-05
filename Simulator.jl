# using DataFrames
# using DataArrays
# using Distributions
#
# include("/Users/austinbean/Desktop/dynhosp/LogitEst.jl")
# include("/Users/austinbean/Desktop/dynhosp/Distance.jl")
# include("/Users/austinbean/Desktop/dynhosp/DemandModel.jl")


#;
#entryprobs = [0.9895, 0.008, 0.0005, 0.002] # [No entry, level1, level2, level3] - not taken from anything, just imposed.
#entrants = [0, 1, 2, 3]
#fields = 7; # if fields updated, update reshaping of state history
#T = 100;

#=
# TESTING --- To run a test, replace the value in the first bracket in mkt_fips = yearins[ ][1], and choose a year.

# Import Data
dataf = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
notmissing = findin(isna(dataf[:fipscode]), false);
dataf = dataf[notmissing, :];
yearins = [ [x; findfirst(dataf[:fipscode], x); findlast(dataf[:fipscode], x ); unique( dataf[findfirst(dataf[:fipscode], x):findlast(dataf[:fipscode], x ) , :year]  ) ] for x in unique(dataf[:fipscode])  ]
mkt_fips = yearins[10][1]
year = 2011
fids = convert(Vector, sort!(unique(dataf[(dataf[:,:fipscode].==mkt_fips)&(dataf[:, :year].==year),:fid])))

# Handle missing dataframe elements::
for i in names(dataf)
    if ( typeof(dataf[i]) == DataArrays.DataArray{Float64,1} )
        dataf[isna(dataf[i]), i] = 0
    elseif (typeof(dataf[i]) == DataArrays.DataArray{Int64,1})
        dataf[isna(dataf[i]), i] = 0
    elseif typeof(dataf[i]) == DataArrays.DataArray{ByteString,1}
        dataf[isna(dataf[i]), i] = "NONE"
  end
    if sum(size(dataf[isna(dataf[i]), i]))>0
    print(i, "\n")
  end
end

data = convert(Matrix, dataf)


# Load the people
people = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);

# Enumerate all of the types::
a = Set()
for el in people.columns
  push!(a, typeof(el))
end

# This is needed to clean out the missing values among fids.  Changes them to 0.
for i in names(people)
  if typeof(people[i]) != DataArrays.DataArray{ByteString,1}
    people[isna(people[i]), i] = 0
  elseif typeof(people[i]) == DataArrays.DataArray{ByteString,1}
    people[isna(people[i]), i] = "NONE"
  end
end

peoples = convert(Matrix, peoplesub)


modcoeffs = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Model.csv", header = true);
distance_c = modcoeffs[1, 2]
distsq_c = modcoeffs[2, 2]
neoint_c = modcoeffs[3, 2]
soloint_c = modcoeffs[4, 2]
closest_c = modcoeffs[5, 2]
distbed_c = modcoeffs[6, 2]

demandmodelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]

# Find the right people:
peoplesub = people[fidfinder(convert(Array, fids)', people, "people"),:]


# To reset for repeated simulations:: (This eliminates entrants, all of which have negative id's)
dataf = dataf[dataf[:id].>= 0, :]

Simulator(dataf, peoplesub, "peoplesub", year, mkt_fips, demandmodelparameters)

To deal with missing values in the hospital dataframe::
for i in names(dataf)
    if ( typeof(dataf[i]) == DataArrays.DataArray{Float64,1} )
        dataf[isna(dataf[i]), i] = 0
    elseif (typeof(dataf[i]) == DataArrays.DataArray{Int64,1})
        dataf[isna(dataf[i]), i] = 0
    elseif typeof(dataf[i]) == DataArrays.DataArray{ByteString,1}
        dataf[isna(dataf[i]), i] = "NONE"
  end
    if sum(size(dataf[isna(dataf[i]), i]))>0
    print(i, "\n")
  end
end


=#


# If I can convert the relevant pieces of the dataframe to an array in the first part, then I won't
# need to hold onto that for the rest of the function - it will probably get much faster.  Something like:
#   mat1 = [ convert(Vector{Float64}, peo[ind[1]]) convert(Vector{Float64}, peo[ind[2]]) convert(Vector{Float64}, peo[ind[3]]) convert(Vector{Float64}, peo[ind[4]]) convert(Vector{Float64}, peo[ind[5]]) convert(Vector{Float64}, peo[ind[6]])]*modelparameters' + rand(d, siz)
# If I can logically index into the array easily then everything should be fine.

# What indices do I use and how do I get them?
# dataf = DataFrames.readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Transition Probabilities.csv", header = true);
# people = DataFrames.readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);

dataf = DataFrames.readtable(pathdata*"TX Transition Probabilities.csv", header = true);
people = DataFrames.readtable(pathpeople*"TX 2005 Individual Choices.csv", header = true);


fipscodeloc = dataf.colindex.lookup[:fipscode]
yearloc = dataf.colindex.lookup[:year]
level1_hospitals0loc = dataf.colindex.lookup[:level1_hospitals0]
level2solo_hospitals0loc = dataf.colindex.lookup[:level2solo_hospitals0]
level3_hospitals0loc = dataf.colindex.lookup[:level3_hospitals0]
fidloc = dataf.colindex.lookup[:fid]
act_intloc = dataf.colindex.lookup[:act_int]
act_sololoc = dataf.colindex.lookup[:act_solo]
lev105loc = dataf.colindex.lookup[:lev105]
lev205loc = dataf.colindex.lookup[:lev205]
lev305loc = dataf.colindex.lookup[:lev305]
lev1515loc = dataf.colindex.lookup[:lev1515]
lev2515loc = dataf.colindex.lookup[:lev2515]
lev3515loc = dataf.colindex.lookup[:lev3515]
lev11525loc = dataf.colindex.lookup[:lev11525]
lev21525loc = dataf.colindex.lookup[:lev21525]
lev31525loc = dataf.colindex.lookup[:lev31525]
v15loc = dataf.colindex.lookup[:v15]
v16loc = dataf.colindex.lookup[:v16]
facilityloc = dataf.colindex.lookup[:facility]
idloc = dataf.colindex.lookup[:id]
locationloc = dataf.colindex.lookup[:location]
cityloc = dataf.colindex.lookup[:city]
firstyearloc = dataf.colindex.lookup[:firstyear]

TotalBeds1loc = people.colindex.lookup[:TotalBeds1]
TotalBeds2loc = people.colindex.lookup[:TotalBeds2]


function Simulator(data::Matrix, peoplesub::Matrix, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2}; T = 100, sim_start = 2, fields = 7, neighbors_start = 108, entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002])
#  if year > 2012
#    return "Years through 2012 only"
#  end
  marketyear = (data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year) # It is a major speed up to do this once.
  level1 = data[marketyear,level1_hospitals0loc][1]
  level2 = data[marketyear,level2solo_hospitals0loc][1]
  level3 = data[marketyear,level3_hospitals0loc][1]
  fids = convert(Vector{Int64}, sort!(unique(data[marketyear,fidloc])))
  state_history = [zeros(1, fields*size(fids, 1)) 1 0 0 0; zeros(T, fields*(size(fids, 1)) + 4)]

  if !( (size(fids)[1])*fields + 4 == size(state_history)[2])
#    println("size of fids: ", size(fids), " size of fields: ", fields, " size of state_history: ", size(state_history))
#    println("first condition: ", (size(fids)[1])*fields + 4, " second condition: ",size(state_history) )
    return "Dims of state_history incorrect"
  end
  # Writes the values to the first row of the state history
  for n in 1:size(fids,1)
    el = fids[n]
    a = ((data[:,fidloc].==el)&marketyear)
    state_history[1, (n-1)*fields + 1] = el # Change
    state_history[1, (n-1)*fields + 2] = data[a,act_intloc][1]
    state_history[1, (n-1)*fields + 3] = data[a,act_sololoc][1]
    state_history[1, (n-1)*fields + 4] = 1 #probability is 1 for the first action
    state_history[1, (n-1)*fields + 5] = 10 # No action at the first period?  Or should it be 10?
    #state_history[1, (n-1)*fields + 6] = # whatever demand is # Not recorded here - recorded below.
    state_history[1, (n-1)*fields + 7] = 0 #record the perturbation 0/1
  end
  # Record aggregte initial values
  state_history[1, (size(fids)[1])*fields+1] = level1 ;
  state_history[1, (size(fids)[1])*fields+2] = level2 ;
  state_history[1, (size(fids)[1])*fields+3] = level3 ;
  state_history[1, (size(fids)[1])*fields+4] = 1; # initial probability.
  # Compute initial demand here:
  peoplesub =rowchange(state_history[1,:], peoplesub)
  #DemandModel(people::DataFrame, frname::ASCIIString, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; maxfid = 11, ent_length = 6 )
  emp_arr = Array{Float64, 2}()
  realized_d = countmap(DemandModel(peoplesub, demandmodelparameters, emp_arr)) # maps chosen hospitals to counts.
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
  end
  # Start simulation here:
  total_entrants = Array{Float64, 2}() # define outside the loop, only for the first round.
  for i = sim_start:T+1
    marketyear = (data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year)
      fids = sort!(unique(data[marketyear,fidloc])) # needs to be updated each round to catch entrants
        for fid in 1:size(fids,1) # this has to be handled separately for each hospital, due to the geography issue
          el = fids[fid] # takes the fid
          a = ((data[:,fidloc].==el)&marketyear)
  #        if sum(a) > 1
  #          println("two entries for ", el, " ", year, " ", mkt_fips)
  #        end
          prev_state = (state_history[i-1, (fid-1)*fields + 2], state_history[i-1, (fid-1)*fields + 3])
  #        print("Part 1 \n")
  #        print("Previous state: ", prev_state, "\n")
          if prev_state ==  (0,0)                     # ((data[a,act_intloc][1], data[a,act_sololoc][1]) == (0,0)) # level 1, actions:
            probs1 = logitest((0,0), level1, level2, level3, [data[a,lev105loc][1]; data[a,lev205loc][1]; data[a,lev305loc][1]; data[a,lev1515loc][1]; data[a,lev2515loc][1]; data[a,lev3515loc][1]; data[a,lev11525loc][1]; data[a,lev21525loc][1]; data[a,lev31525loc][1]] )
            if probs1 == ValueException
    #          println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
              break
            end
            # Draw action:
            action1 = sample([10, 2, 1, 11] ,WeightVec([probs1[1], probs1[2], probs1[3], probs1[4]]))
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
            end
            # Set own distance counts to 0 for all categories
            (data[a,lev105loc], data[a,lev205loc], data[a,lev305loc], data[a,lev1515loc], data[a,lev2515loc], data[a,lev3515loc], data[a,lev11525loc], data[a,lev21525loc], data[a,lev31525loc]) = zeros(1,9)
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
            if probs2 == ValueException
              println("Value Exception at ", level1, " ", level2, " ", level3, " ", mkt_fips, " year", year)
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
            if probs3 == ValueException
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
        end
        # Measure distances to neighbors
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
      entrantsp = WeightVec(entryprobs)
      newentrant = sample(entrants, entrantsp)
      entrantout = [newentrant, entrypairs[findfirst(entrypairs[:,1], newentrant), 2]]'
      b = ((data[ :,fipscodeloc].== mkt_fips)&(data[ :,yearloc].==year)); # Market-year observations, whether exited or not.
        if newentrant> 0
    #      println("Entry occurred")
          # Eventually fix the fact that the neighbors here are not going to be exactly right.
          # I need to add the fact that the entrants are not being recorded in the "neighbors" section so no one is counting distances to them.
          ent_lat = mean((data[b,v15loc])) + rand(Normal(0, 0.1), 1) # 0.1 degrees latitude should be about 6-7 miles.
          ent_lon = mean((data[b,v16loc])) + rand(Normal(0, 0.1), 1)
          newrow = data[b,:][1,:] # create new row, duplicating existing.  Takes first row of current
          newrow[facilityloc] = convert(UTF8String, "Entrant $mkt_fips $year")
          newrow[fidloc] = sample(collect((maximum(data[:,fidloc])+1):(maximum(data[:,fidloc])+5))) # new facility will always have largest fid
          newrow[idloc] = -sample(collect((maximum(data[:,idloc])+1):(maximum(data[:,idloc])+5))) #
          newrow[fipscodeloc] = mkt_fips
          newrow[locationloc] = convert(UTF8String, "entrant - see v15 v16")
          newrow[cityloc] = convert(UTF8String, "Entrant - unspecified")
          newrow[firstyearloc] = year
          newrow[v15loc] = ent_lat
          newrow[v16loc] = ent_lon
          # Take the size as the mean bed number from neighboring hospitals.  There is no field for this in dataf, unfortunately.
          # Just compute the mean over all beds in the state?  This needs to be fixed later.
          entrantbeds = convert(Int, floor(mean( unique(vcat(unique(peoplesub[ peoplesub[TotalBeds1loc].>0 ,TotalBeds1loc]), unique(peoplesub[ peoplesub[TotalBeds2loc].>0 ,TotalBeds2loc])) ))) )
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
              if valfn > 0
                nbloc = neighbors_start + 2*valfn
                data[row, nbloc ]= newrow[fidloc]
                data[row, nbloc+1]= distance(ent_lat[1], ent_lon[1], newrow[v15loc][1], newrow[v16loc][1])
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
          data = [data; newrow]
          # append value to fids
          push!(fids, newrow[fidloc][1])
          # Reshape state history: fid, solo state, int state, probability of choice,  action chosen, XXXX demand, perturbed. [newrow[:fid], 999, 999, 0, 1, 0]
          #    println("reshaping state history")
          # The problem is the size computation right here - figure it out.  I'm not resizing this right.  The entry condition is not correct now.
          # OLD:  vcat(hcat(state_history[1:i,1:end-4], repmat([newrow[fidloc][1] 999 999 1 0 0 0], i, 1), state_history[1:i, end-3:end]), zeros((T-i+1), size(fids,1)*fields+4 ))
          # Want: vcat( hcat(state_history[1:i-1, 1:end-4], repmat([newrow[fidloc][1] 999 999 1 0 0 0], i-1, 1), state_history[1:i-1,end-3:end]), hcat(state_history[i,1:end-4], [newrow[fidloc] newrow[act_sololoc] newrow[act_intloc] entrantout[2] 0 0 0], state_history[i, end-3:end]) ,zeros((T-i+1), (size(fids,1)+1)*fields+4) ))
          state_history = vcat( hcat(state_history[1:i-1, 1:end-4], repmat([newrow[fidloc][1] 999 999 1 0 0 0], i-1, 1), state_history[1:i-1,end-3:end]), hcat(state_history[i,1:end-4], [newrow[fidloc] newrow[act_sololoc] newrow[act_intloc] entrantout[2] 0 0 0], state_history[i, end-3:end]) ,zeros((T-i+1), size(fids,1)*fields+4) )
        end
        # Aggregate Probability of Action:
      tprob = 1
#      print("Part 4", "\n")
      for els in 4:fields:(size(state_history[i,:], 2)-4) # 4 is the relevant column
        # if (state_history[i,els-2] == -999) | (state_history[i,els-1] == -999)
        #   print("should be the fid ", state_history[i,els-3], "\n")
        #   print("And the row: ", i, "\n")
        #   print("And the elem: ", els, "\n")
        #   print("The previous state row: ", i-1,"  ", showall(state_history[i-1,:]'), "\n")
        #   print("******", "\n")
        #   print("The current state row: ", i, "   ", showall(state_history[i,:]'), "\n")
        #   print("\n")
        # end
        if (state_history[i,els] <= 0) | (state_history[i,els]>1)
            # print("elements < 0? ", state_history[i,els] <= 0, "\n")
            # print("the element: ", state_history[i,els], " and the i ", i, "\n")
            # print("nearby elements: ", state_history[i, els-3:els+3], "\n")
            # print("elements > 1? ", state_history[i,els]>1, "\n")
    #        print("Bad probability at row ", i, " ", els, " ", state_history[i,els-3 ], " prob ", state_history[i,els], "\n" )
            # print(showall(state_history[i-1,:]'), "\n")
            # print("******", "\n")
            # print(showall(state_history[i,:]'), "\n")
            # return state_history
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
      # Since people sub is running twice, it's probably taking 10 - 12 seconds per iteration.  Must be sped up.
      peoplesub = rowchange(state_history[i-1,:], peoplesub)
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
  return state_history
end
