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

# Duplicates the above but takes fidvect::Array{Int64, 1} if necessary.

function fidfinder(fidvect::Array{Int64, 1}, choices::DataFrame, frname::ASCIIString; maxfid = 11)
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


function rowchange(staterow::Array{Float64,2}, choicerow::DataFrame; endfields_state = 4, fields_state = 7, fields_people = 16, endfields_people = 7)
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
    peoplefids =  unique([choicerow[x][1] for x in 2:fields_people:size(choicerow)[2]-(endfields_people) ]) # collects all fids in the person's choice set
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
        elseif (staterow[findfirst(staterow, el)+1] == -999) | (staterow[findfirst(staterow, el)+2] == -999)
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

function EntrantsU(peo::DataFrame, entrants::Array{Float64, 2}, modelparameters::Array{Float64, 2}; persloc = [183 184], entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize))
  siz = size(peo,1)
  rands = rand(d, siz, entnum)
  entvals = zeros(siz, entnum)
  plocs = [convert(Vector{Float64}, peo[persloc[1]]) convert(Vector{Float64}, peo[persloc[2]])] # this format is: LATITUDE, LONGITUDE
  for j = 1:entnum
    for i = 1:siz # ENTRANT FORMAT ends with [Latitude, Longitude]
      d1 = distance(entrants[6*j-1], entrants[6*j], plocs[i,1], plocs[i,2])
      if d1 < 50
        entvals[i,j] += d1*modelparameters[1]
        entvals[i,j] += ((d1)^2)*modelparameters[2]
        entvals[i,j] += entrants[6*j-3]*modelparameters[3]
        entvals[i,j] += entrants[6*j-4]*modelparameters[4]
        entvals[i,j] += 0*modelparameters[5] # this is specifying that "closest" is always 0 for entrants.  It can be fixed, but would be really annoying.
        entvals[i,j] += ((d1)*entrants[6*j-2]/100)*modelparameters[6]
      else
        entvals[i,j] = -999 # set value to large negative number when distance is too large: won't be chosen
      end
    end
  end
  return maximum(entvals + rands, 2) # note that due to the randomization, this will generally not return -999, but -999 + rand
end



function DemandModel(peo::DataFrame, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), persloc = [183 184] , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162] )
# constants/outputs/setup
  outp = zeros(siz)
  entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
# Computed utilities + error
  mat1 = [ convert(Vector{Float64}, peo[ind[1]]) convert(Vector{Float64}, peo[ind[2]]) convert(Vector{Float64}, peo[ind[3]]) convert(Vector{Float64}, peo[ind[4]]) convert(Vector{Float64}, peo[ind[5]]) convert(Vector{Float64}, peo[ind[6]])]*modelparameters' + rand(d, siz)
  mat2 = [ convert(Vector{Float64}, peo[iind[1]]) convert(Vector{Float64}, peo[iind[2]]) convert(Vector{Float64}, peo[iind[3]]) convert(Vector{Float64}, peo[iind[4]]) convert(Vector{Float64}, peo[iind[5]]) convert(Vector{Float64}, peo[iind[6]])]*modelparameters' + rand(d, siz)
  mat3 = [ convert(Vector{Float64}, peo[iiind[1]]) convert(Vector{Float64}, peo[iiind[2]]) convert(Vector{Float64}, peo[iiind[3]]) convert(Vector{Float64}, peo[iiind[4]]) convert(Vector{Float64}, peo[iiind[5]]) convert(Vector{Float64}, peo[iiind[6]])]*modelparameters' + rand(d, siz)
  mat4 = [ convert(Vector{Float64}, peo[ivnd[1]]) convert(Vector{Float64}, peo[ivnd[2]]) convert(Vector{Float64}, peo[ivnd[3]]) convert(Vector{Float64}, peo[ivnd[4]]) convert(Vector{Float64}, peo[ivnd[5]]) convert(Vector{Float64}, peo[ivnd[6]])]*modelparameters' + rand(d, siz)
  mat5 = [ convert(Vector{Float64}, peo[vnd[1]]) convert(Vector{Float64}, peo[vnd[2]]) convert(Vector{Float64}, peo[vnd[3]]) convert(Vector{Float64}, peo[vnd[4]]) convert(Vector{Float64}, peo[vnd[5]]) convert(Vector{Float64}, peo[vnd[6]])]*modelparameters' + rand(d, siz)
  mat6 = [ convert(Vector{Float64}, peo[vind[1]]) convert(Vector{Float64}, peo[vind[2]]) convert(Vector{Float64}, peo[vind[3]]) convert(Vector{Float64}, peo[vind[4]]) convert(Vector{Float64}, peo[vind[5]]) convert(Vector{Float64}, peo[vind[6]])]*modelparameters' + rand(d, siz)
  mat7 = [ convert(Vector{Float64}, peo[viind[1]]) convert(Vector{Float64}, peo[viind[2]]) convert(Vector{Float64}, peo[viind[3]]) convert(Vector{Float64}, peo[viind[4]]) convert(Vector{Float64}, peo[viind[5]]) convert(Vector{Float64}, peo[viind[6]])]*modelparameters' + rand(d, siz)
  mat8 = [ convert(Vector{Float64}, peo[viiind[1]]) convert(Vector{Float64}, peo[viiind[2]]) convert(Vector{Float64}, peo[viiind[3]]) convert(Vector{Float64}, peo[viiind[4]]) convert(Vector{Float64}, peo[viiind[5]]) convert(Vector{Float64}, peo[viiind[6]])]*modelparameters' + rand(d, siz)
  mat9 = [ convert(Vector{Float64}, peo[ixnd[1]]) convert(Vector{Float64}, peo[ixnd[2]]) convert(Vector{Float64}, peo[ixnd[3]]) convert(Vector{Float64}, peo[ixnd[4]]) convert(Vector{Float64}, peo[ixnd[5]]) convert(Vector{Float64}, peo[ixnd[6]])]*modelparameters' + rand(d, siz)
  mat10 = [ convert(Vector{Float64}, peo[xnd[1]]) convert(Vector{Float64}, peo[xnd[2]]) convert(Vector{Float64}, peo[xnd[3]]) convert(Vector{Float64}, peo[xnd[4]]) convert(Vector{Float64}, peo[xnd[5]]) convert(Vector{Float64}, peo[xnd[6]])]*modelparameters' + rand(d, siz)
  mat11 = [ convert(Vector{Float64}, peo[xind[1]]) convert(Vector{Float64}, peo[xind[2]]) convert(Vector{Float64}, peo[xind[3]]) convert(Vector{Float64}, peo[xind[4]]) convert(Vector{Float64}, peo[xind[5]]) convert(Vector{Float64}, peo[xind[6]])]*modelparameters' + rand(d, siz)
  if size(entrants, 2) > 1
    allfids = [convert(Vector{Int64}, peo[fidnd[1]]) convert(Vector{Int64}, peo[fidnd[2]]) convert(Vector{Int64}, peo[fidnd[3]]) convert(Vector{Int64}, peo[fidnd[4]]) convert(Vector{Int64}, peo[fidnd[5]]) convert(Vector{Int64}, peo[fidnd[6]]) convert(Vector{Int64}, peo[fidnd[7]]) convert(Vector{Int64}, peo[fidnd[8]]) convert(Vector{Int64}, peo[fidnd[9]]) convert(Vector{Int64}, peo[fidnd[10]]) convert(Vector{Int64}, peo[fidnd[11]]) repmat(entfids, siz, 1)]
    entval = EntrantsU(people, entrants, modelparameters)
    for i = 1:size(peo, 1)
      outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i] entval[i]])]
    end
  else # there are no entrants, so "entrants" above is [0.0]', which has size(entrants, 2) == 1
    allfids = [convert(Vector{Int64}, peo[fidnd[1]]) convert(Vector{Int64}, peo[fidnd[2]]) convert(Vector{Int64}, peo[fidnd[3]]) convert(Vector{Int64}, peo[fidnd[4]]) convert(Vector{Int64}, peo[fidnd[5]]) convert(Vector{Int64}, peo[fidnd[6]]) convert(Vector{Int64}, peo[fidnd[7]]) convert(Vector{Int64}, peo[fidnd[8]]) convert(Vector{Int64}, peo[fidnd[9]]) convert(Vector{Int64}, peo[fidnd[10]]) convert(Vector{Int64}, peo[fidnd[11]])]
    for i = 1:size(peo, 1)
      outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i]])]
    end
  end
  return outp
end


# use countmap(choicemade) to count the results (!)  So easy.



#allents = zeros(size(people)[1], convert(Int, maxfid + floor((size(entrants)[2])/6)))

function DemandMod(people::DataFrame, frname::ASCIIString, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; maxfid = 11, ent_length = 6 )
# Performance: 19.633987 seconds (3.01 M allocations: 233.659 MB, 1.04% gc time) (using whole 2005 Q1 dataset)
  allents = zeros(size(people)[1], maxfid+ convert(Int, floor(size(entrants)[2]/ent_length)))
  for i in 1:maxfid
    shock = rand(d, size(people)[1])
    expr = parse("$frname[:distance$i].data*distance_c + $frname[:distsq$i].data*distsq_c + $frname[:SoloIntermediate$i].data*soloint_c + $frname[:NeoIntensive$i].data*neoint_c + $frname[:closest$i].data*closest_c + $frname[:dist_bed$i].data*distbed_c + shock  ")
    allents[:,i] =  eval(expr)
  end
  votes = zeros(size(people)[1], 1);
  for row in 1:size(allents)[1]
    votes[row] = eval(parse("people[$row,:fid"*string(indmax(allents[i,:]))*"]"))
  end
  return countmap(votes)
end


#Performance: 0.510989 seconds (5.64 M allocations: 213.967 MB, 68.09% gc time)






















#end
