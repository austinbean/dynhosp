# Use the parameters of the demand model to estimate choices:

# Version: 05 09 16

using DataFrames
using DataArrays
using Distributions

# Now in use in Main.jl

function fidfinder(fidvect::Array{Int64, 2}, choices::DataFrame, frname::ASCIIString; maxfid = 11)
    #=
      This function generates a vector of Booleans which index the rows in the
      dataframe in which the individual has some hospital with fid in fidvect
      as an option.  It starts with all falses and iteratively takes subset | (result)
      which will be true when the result expression is true.  It operates on a whole DataFrame
      The function takes as one argument the name of the frame "choices" (frname) as a string.
    =#
    subset = falses(size(choices)[1])
      for k in 1:size(fidvect)[1]
        targ = fidvect[k]
        for j = 1:maxfid
           subex = parse(frname*"[:fid$j].==$targ")
           true_v = eval(subex)
      #     print( size(true_v), "   ") # for testing purposes
           subset = subset | true_v
        end
      end
    return subset
end

# The next function will find values in the one-row dataframe element given a list of symbols


function rowfindfid(targ::DataFrame, value::Int64; vals = [:fid1, :fid2, :fid3, :fid4, :fid5, :fid6, :fid7, :fid8, :fid9, :fid10, :fid11] )
  #=
    This function takes a single row of a dataframe "targ" and finds the value in "value"
    by looking in the columns in "vals".  Those items in vals must be a list of symbols, so
    [ :fid1, :fid2, ..., :fidN ] (commas necessary).  If the DataFrame is too big, it Returns
    the RowSizeError exception defined above.  To search specific fids, must call function as
    rowfindfid(targ, value, vals = [:fidx, :fidy])
  =#
  if size(targ)[1] > 1
    return RowSizeError #defined in Main.jl
  end
  numb = 0
  index = :ident
  for i in vals # what about NAs?  All of the fid columns are zeros.
      if targ[i][1] == value # the second [1] is needed to access the *value* of the element, rather than the 1 element DataFrame
          numb = targ.colindex.lookup[i] # returns the *column number* of the element i
          index =  i
      end
  end
  return numb, index
end


function rowchange(staterow::Array{Float64,2}, choicerow::DataFrame; endfields_state = 4, fields_state = 7, fields_people = 15, endfields_people = 7)
  #=
     This function should take a row of the state history (staterow), and a row of
     the choices (choicerow) and:
     1.  determines the number of fids in the staterow ✓
     2.  Determines the number of fids in the choicerow ✓
     3.  When a fid in the staterow matches a fid in the choicerow, map the values
         from the staterow to the choicerow
   Notes - need to do something special for entrants.
   Can check if fid sets are overlapping - change those fids which are
   Once I know this, I also need to check whether the new hospital is the closest.

   staterow has the form: [ fid, act_solo, act_int, choice prob, action taken, demand realized, perturbed] × (# facilities)  ⋃ [ level1's total, level2's total, level3's total, aggregate prob]
   choicerow has the form: [identity, fid, facility, NeoIntensive, TotalDeliveries, Transfers Out No NICU, Transfers In Has NICU, Not For Profit Status (#), Solo Intermediate, distance, Is Closest?, Selected?, NFP ?, distance × bed, distance²] × (# facilities) ⋃ [Patient Zip, CMS MDC, APR MDC, CMS DRG, APR DRG, Zip Lat, Zip Long]

  =#
    # Collects the fids which are in the market
    mktnumfids = unique(((size(staterow)[2])-endfields_state)/fields_state ) # number of facilities
    mktfids = [ el for el in staterow[1,1:fields_state:end-endfields_state]] # Collects the fids in the market
    # Collects the fids which are in the choice set
    peoplefids =  unique([choicerow[x][1] for x in 2:fields_people:size(choicerow)[2]-endfields_people ]) # collects all fids in the person's choice set
    peoplenumfids = unique(sum(peoplefids.>0)) # Counts the number of unique facilities (fid > 0) in the choice set (missing facilities have fid = 0, rather than NA)

    # Takes the values of market fids which are in the choice row (only these must be changed)
    change_fids = intersect(peoplefids, mktfids)
    if sum(size(change_fids))== 0
      return choicerow
    else # intersection is nonzero
      for el in change_fids
        # Here - findfirst(staterow, el)
        # If that + 1 == -999 and that + 2 = -999
        # Set that fid to 0 (treats the hospital as missing)
        # Then need to reload "people" later.
        if staterow[findfirst(staterow, el)+1] != -999
          el = convert(Int64, el)
          (loc, symb) = rowfindfid(choicerow, el) #finds the fid in the row, or returns 0 if absent
          if loc != 0
            fid_num = replace(string(symb), r"[fid]", "") # takes the name of the symbol (:fid#), converts to string, "fid#", removes 'fid', obtains "#" as string.
            # Change NeoIntensive in the choice row to the value in the state row
            neo = Symbol("NeoIntensive"*fid_num)
            choicerow[neo] = convert(Int64, staterow[findfirst(staterow, el)+2])
            # Change SoloIntermeidate in the choice row to the value in the state row
            solo = Symbol("SoloIntermediate"*fid_num)
            choicerow[solo] = convert(Int64, staterow[findfirst(staterow, el)+1])
          else
            print(el, " not found in row ")
          end
        elseif staterow[findfirst(staterow, el)+1] == -999
          el = convert(Int64, el)
          (loc, symb) = rowfindfid(choicerow, el) #finds the fid in the row, or returns 0 if absent
          if loc != 0
            choicerow[symb] = 0 # reassign the value of fid to be zero so that demand cannot be computed for an exited hospital
          end
        end
      end
    end
    return choicerow
end


dist_μ = 0;
dist_σ = 1;
dist_ξ = 0;
srand(123)
d = GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# I do need the constant:
γ = eulergamma;





function DemandModel(people::DataFrame, frname::ASCIIString, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; maxfid = 11, ent_length = 6 )
  #=
    The goal for this function is to -
      - take the whole set of people, compute the deterministic components of utility, add the random shock, find the maximizer
      - count the number maximized by fid: this will be the demand.
      Performance: .518402 seconds (32.05 M allocations: 604.392 MB, 0.76% gc time)
      The first two arguments are a dataframe (people) and the NAME of that dataframe (frname) as a string.  This is important at the end
      - Entrants - measure length, computes number of entrants.
      - Format of entrants is: [fid act_int act_solo entrantbeds ent_lat ent_lon] x (# entrants)
      - Remember: distbed is distance * beds/100
  =#
  choicemade = zeros(Int64, size(people)[1], 1)
  distance_c  = modelparameters[1]
  distsq_c  = modelparameters[2]
  neoint_c  = modelparameters[3]
  soloint_c  = modelparameters[4]
  closest_c  = modelparameters[5]
  distbed_c = modelparameters[6]
  for i in 1:size(people)[1]
    shock = rand(d, maxfid)
    if people[i,:fid1] != 0
      val1 = people[i,:distance1]*distance_c + people[i,:distsq1]*distsq_c + people[i,:SoloIntermediate1]*soloint_c + people[i,:NeoIntensive1]*neoint_c + people[i,:closest1]*closest_c + people[i,:dist_bed1]*distbed_c + shock[1]
      dist1 = people[i,:distance1]
    else
      val1 = -10^10
      dist1 = 999
    end
    if people[i,:fid2] != 0
      val2 = people[i,:distance2]*distance_c + people[i,:distsq2]*distsq_c + people[i,:SoloIntermediate2]*soloint_c + people[i,:NeoIntensive2]*neoint_c + people[i,:closest2]*closest_c + people[i,:dist_bed2]*distbed_c + shock[2]
      dist2 = people[i,:distance2]
    else
      val2 = -10^10
      dist2 = 999
    end
    if people[i,:fid3] != 0
      val3 = people[i,:distance3]*distance_c + people[i,:distsq3]*distsq_c + people[i,:SoloIntermediate3]*soloint_c + people[i,:NeoIntensive3]*neoint_c + people[i,:closest3]*closest_c + people[i,:dist_bed3]*distbed_c + shock[3]
      dist3 = people[i,:distance3]
    else
      val3 = -10^10
      dist3 = 999
    end
    if people[i,:fid4] != 0
      val4 = people[i,:distance4]*distance_c + people[i,:distsq4]*distsq_c + people[i,:SoloIntermediate4]*soloint_c + people[i,:NeoIntensive4]*neoint_c + people[i,:closest4]*closest_c + people[i,:dist_bed4]*distbed_c + shock[4]
      dist4 = people[i,:distance4]
    else
      val4 = -10^10
      dist4 = 999
    end
    if people[i,:fid5] != 0
      val5 = people[i,:distance5]*distance_c + people[i,:distsq5]*distsq_c + people[i,:SoloIntermediate5]*soloint_c + people[i,:NeoIntensive5]*neoint_c + people[i,:closest5]*closest_c + people[i,:dist_bed5]*distbed_c + shock[5]
      dist5 = people[i,:distance5]
    else
      val5 = -10^10
      dist5 = 999
    end
    if people[i,:fid6] != 0
      val6 = people[i,:distance6]*distance_c + people[i,:distsq6]*distsq_c + people[i,:SoloIntermediate6]*soloint_c + people[i,:NeoIntensive6]*neoint_c + people[i,:closest6]*closest_c + people[i,:dist_bed6]*distbed_c + shock[6]
      dist6 = people[i,:distance6]
    else
      val6 = -10^10
      dist6 = 999
    end
    if people[i,:fid7] != 0
      val7 = people[i,:distance7]*distance_c + people[i,:distsq7]*distsq_c + people[i,:SoloIntermediate7]*soloint_c + people[i,:NeoIntensive7]*neoint_c + people[i,:closest7]*closest_c + people[i,:dist_bed7]*distbed_c + shock[7]
      dist7 = people[i,:distance7]
    else
      val7 = -10^10
      dist7 = 999
    end
    if people[i,:fid8] != 0
      val8 = people[i,:distance8]*distance_c + people[i,:distsq8]*distsq_c + people[i,:SoloIntermediate8]*soloint_c + people[i,:NeoIntensive8]*neoint_c + people[i,:closest8]*closest_c + people[i,:dist_bed8]*distbed_c + shock[8]
      dist8 = people[i,:distance8]
    else
      val8 = -10^10
      dist8 = 999
    end
    if people[i,:fid9] != 0
      val9 = people[i,:distance9]*distance_c + people[i,:distsq9]*distsq_c + people[i,:SoloIntermediate9]*soloint_c + people[i,:NeoIntensive9]*neoint_c + people[i,:closest9]*closest_c + people[i,:dist_bed9]*distbed_c + shock[9]
      dist9 = people[i,:distance9]
    else
      val9 = -10^10
      dist9 = 999
    end
    if people[i,:fid10] != 0
      val10 = people[i,:distance10]*distance_c + people[i,:distsq10]*distsq_c + people[i,:SoloIntermediate10]*soloint_c + people[i,:NeoIntensive10]*neoint_c + people[i,:closest10]*closest_c + people[i,:dist_bed10]*distbed_c + shock[10]
      dist10 = people[i,:distance10]
    else
      val10 = -10^10
      dist10 = 999
    end
    if people[i,:fid11] != 0
      val11 = people[i,:distance11]*distance_c + people[i,:distsq11]*distsq_c + people[i,:SoloIntermediate11]*soloint_c + people[i,:NeoIntensive11]*neoint_c + people[i,:closest11]*closest_c + people[i,:dist_bed11]*distbed_c + shock[11]
      dist11 = people[i,:distance11]
    else
      val11 = -10^10
      dist11 = 999
    end
    choice = 0;
    if maximum(size(entrants)) <= 1 # there isn't a two-d array with size 0.
      chosen = indmax([val1 val2 val3 val4 val5 val6 val7 val8 val9 val10 val11]) # returns the *index* of the max element in the collection
    else # at least one entrant
      cval = maximum([val1 val2 val3 val4 val5 val6 val7 val8 val9 val10 val11])
      mindist = minimum([dist1 dist2 dist3 dist4 dist5 dist6 dist7 dist8 dist9 dist10 dist11])
      chosen = indmax([val1 val2 val3 val4 val5 val6 val7 val8 val9 val10 val11])
      entfids = [entrants[x] for x in 1:ent_length:maximum(size(entrants))]
      for k = 1:size(entfids)[1]
        efid = entfids[k]
        eind = findfirst(entrants, efid)
        # generate the distance x beds interaction
        dist = distance(people[i,:ZIP_LAT], people[i, :ZIP_LONG], entrants[eind+4], entrants[eind+5])
        if dist < 50
          if dist < mindist
            cl_ind = 1
          else
            cl_ind = 0
          end
          val = distance_c*dist + distsq_c*dist^2 + neoint_c*entrants[eind+1] + soloint_c*entrants[eind+2] + closest_c*cl_ind + distbed_c*dist*entrants[eind+3]/100 + rand(d, 1) # need to fix closest later
          # Strategy - compare to chosen, if larger, take that.  Return chosen = 12 later can do this sequentially.
          if val[1] > cval[1]
            chosen = 99
            choice = k
            cval = val # reassign
          end
        end
      end
    end
    if chosen <= 11
      choicemade[i] = convert(Int64, eval(parse(frname*"[$i, :fid$chosen]")))
    else
      choicemade[i] = convert(Int64, entfids[choice])
    end

  end
  return choicemade
end




# use countmap(choicemade) to count the results (!)  So easy.


#=
for i in 1:maxfid
  # This is vectorized - probably slow.  Why not just go by row and see if it's faster?
  expr = parse("vals$i = hcat(people[:fid$i] , people[:distance$i].data*distance_c + people[:distsq$i].data*distsq_c + people[:SoloIntermediate$i].data*soloint_c + people[:NeoIntensive$i].data*neoint_c + people[:closest$i].data*closest_c + people[:dist_bed$i].data*distbed_c  )") # evaluate a tuple
  eval(expr)
end

Performance: 0.510989 seconds (5.64 M allocations: 213.967 MB, 68.09% gc time)
=#





















#end
