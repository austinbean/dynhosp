# Use the parameters of the demand model to estimate choices:

# Version: 05 09 16
#
# using DataFrames
# using DataArrays
# using Distributions

# Now in use in Main.jl

# Indices:
#people = DataFrames.readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);
people = DataFrames.readtable(pathpeople*"TX 2005 1 Individual Choices.csv", header = true);

fid1loc = people.colindex.lookup[:fid1]
fid2loc = people.colindex.lookup[:fid2]
fid3loc = people.colindex.lookup[:fid3]
fid4loc = people.colindex.lookup[:fid4]
fid5loc = people.colindex.lookup[:fid5]
fid6loc = people.colindex.lookup[:fid6]
fid7loc = people.colindex.lookup[:fid7]
fid8loc = people.colindex.lookup[:fid8]
fid9loc = people.colindex.lookup[:fid9]
fid10loc = people.colindex.lookup[:fid10]
fid11loc = people.colindex.lookup[:fid11]

function fidfinder(fidvect::Array{Int64, 2}, choices::Matrix; maxfid = 11)
    #=
      This function generates a vector of Booleans which index the rows in the
      dataframe in which the individual has some hospital with fid in fidvect
      as an option.  It starts with all falses and iteratively takes subset | (result)
      which will be true when the result expression is true.  It operates on a whole DataFrame
      The function takes as one argument the name of the frame "choices" (frname) as a string.
    =#
    subset = falses(size(choices)[1])
      for k in 1:maximum(size(fidvect))
        targ = fidvect[k]
        for j in [fid1loc fid2loc fid3loc fid4loc fid5loc fid6loc fid7loc fid8loc fid9loc fid10loc fid11loc]
           subex = (choices[:,j].==targ)
      #     true_v = eval(subex)
      #     print( size(true_v), "   ") # for testing purposes
           subset = subset | subex
        end
      end
    return subset
end

# Duplicates the above but takes fidvect::Array{Int64, 1} if necessary.
# Timing: @time fidfinder(convert(Array, fids)', people, "people")
#  0.010568 seconds (837 allocations: 4.533 MB)

function fidfinder(fidvect::Array{Int64, 1}, choices::Matrix; maxfid = 11)
    #=
      This function generates a vector of Booleans which index the rows in the
      dataframe in which the individual has some hospital with fid in fidvect
      as an option.  It starts with all falses and iteratively takes subset | (result)
      which will be true when the result expression is true.  It operates on a whole DataFrame
      The function takes as one argument the name of the frame "choices" (frname) as a string.
    =#
      subset = falses(size(choices)[1])
      for k in 1:maximum(size(fidvect))
        targ = fidvect[k]
        for j in [fid1loc fid2loc fid3loc fid4loc fid5loc fid6loc fid7loc fid8loc fid9loc fid10loc fid11loc]
           subex = (choices[:,j].==targ)
           subset = subset | subex
        end
      end
    return subset
end


# The next function will find values in the one-row dataframe element given a list of symbols
# For use in the demand model::
type  RowSizeError  <: Exception end



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

# Large number of allocations: 0.699547 seconds (4.58 M allocations: 167.552 MB)
# This for a matrix with 89000 rows.
function rowchange(staterow::Array{Float64,2}, choicerow::Matrix; choiceintloc = 3, choicesololoc = 8, endfields_state = 4, fields_state = 7, fields_people = 16, endfields_people = 7)
  #=  The first argument is the state history, the second the individual record
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
    mktfids = [ el for el in staterow[1,1:fields_state:end-endfields_state]] # Collects the fids in the market
    # Collects the fids which are in the choice set
    for i in 1:size(choicerow, 1)
      peoplefids =  convert(Array{Int64}, unique([choicerow[i,x][1] for x in 2:fields_people:size(choicerow, 2)-(endfields_people) ])) # collects all fids in the person's choice set
      change_fids = intersect(peoplefids, mktfids) # Takes the values of market fids which are in the choice row (only these must be changed)
      if sum(size(change_fids))!= 0
        for el in change_fids
          fidloc = findfirst(choicerow[i,:], el) # finds the fid in the row of the person's choices
          statefidloc = findfirst(staterow, el) # finds the fid in the state history row
          if (staterow[statefidloc+1] != -999) | (staterow[statefidloc+2] != -999) #check whether firm exited
            choicerow[i,fidloc+choicesololoc] = staterow[statefidloc+1]
            choicerow[i,fidloc+choiceintloc] = staterow[statefidloc+2]
          elseif (staterow[statefidloc+1] == -999) | (staterow[statefidloc+2] == -999)
            choicerow[i,fidloc] = 0
            choicerow[i,fidloc+choicesololoc] = -999
            choicerow[i,fidloc+choiceintloc] = -999
          end
        end
      end
    end
  return choicerow
end

#
# dist_μ = 0;
# dist_σ = 1;
# dist_ξ = 0;
# srand(123)
# d = GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# # I do need the constant:
γ = eulergamma;

# sample entrants1 = [99999 1 0 120 32.96  -96.8385] [newrow[fidloc] newrow[act_intloc] newrow[act_sololoc] entrantbeds ent_lat ent_lon]
# sample entrants2 = [99999 1 0 120 32.96  -96.8385 888888 0 1 120 32.96  -96.8385]

function EntrantsU(peo::Matrix, entrants::Array{Float64, 2}, modelparameters::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, persloc = [183 184], entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize))
  d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
  siz = size(peo,1)
  rands = rand(d, siz, entnum)
  entvals = zeros(siz, entnum)
  plocs = [convert(Vector{Float64}, peo[:,persloc[1]]) convert(Vector{Float64}, peo[:,persloc[2]])] # this format is: LATITUDE, LONGITUDE
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
  utils, inds = findmax(entvals + rands, 2)
  return [utils ind2sub(size(entvals), vec(inds))[2]]  #maximum(entvals + rands, 2) # note that due to the randomization, this will generally not return -999, but -999 + rand
end


# 1.599738 seconds (9.47 M allocations: 355.385 MB, 56.89% gc time)
# 0.533586 seconds (7.71 M allocations: 310.720 MB)
# entrants = Array{Float64, 2}() # for a fast test with no entrants
# That's more costly than I would like
# Experiment with eliminating the conversions, just form the matrices.
function DemandModel(peo::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), persloc = [183 184] , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162] )
  d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# constants/outputs/setup
  outp = zeros(siz)
  entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
# Computed utilities + error
  mat1 = [ convert(Vector{Float64}, peo[:,ind[1]]) convert(Vector{Float64}, peo[:,ind[2]]) convert(Vector{Float64}, peo[:,ind[3]]) convert(Vector{Float64}, peo[:,ind[4]]) convert(Vector{Float64}, peo[:,ind[5]]) convert(Vector{Float64}, peo[:,ind[6]])]*modelparameters' + rand(d, siz)
  mat2 = [ convert(Vector{Float64}, peo[:,iind[1]]) convert(Vector{Float64}, peo[:,iind[2]]) convert(Vector{Float64}, peo[:,iind[3]]) convert(Vector{Float64}, peo[:,iind[4]]) convert(Vector{Float64}, peo[:,iind[5]]) convert(Vector{Float64}, peo[:,iind[6]])]*modelparameters' + rand(d, siz)
  mat3 = [ convert(Vector{Float64}, peo[:,iiind[1]]) convert(Vector{Float64}, peo[:,iiind[2]]) convert(Vector{Float64}, peo[:,iiind[3]]) convert(Vector{Float64}, peo[:,iiind[4]]) convert(Vector{Float64}, peo[:,iiind[5]]) convert(Vector{Float64}, peo[:,iiind[6]])]*modelparameters' + rand(d, siz)
  mat4 = [ convert(Vector{Float64}, peo[:,ivnd[1]]) convert(Vector{Float64}, peo[:,ivnd[2]]) convert(Vector{Float64}, peo[:,ivnd[3]]) convert(Vector{Float64}, peo[:,ivnd[4]]) convert(Vector{Float64}, peo[:,ivnd[5]]) convert(Vector{Float64}, peo[:,ivnd[6]])]*modelparameters' + rand(d, siz)
  mat5 = [ convert(Vector{Float64}, peo[:,vnd[1]]) convert(Vector{Float64}, peo[:,vnd[2]]) convert(Vector{Float64}, peo[:,vnd[3]]) convert(Vector{Float64}, peo[:,vnd[4]]) convert(Vector{Float64}, peo[:,vnd[5]]) convert(Vector{Float64}, peo[:,vnd[6]])]*modelparameters' + rand(d, siz)
  mat6 = [ convert(Vector{Float64}, peo[:,vind[1]]) convert(Vector{Float64}, peo[:,vind[2]]) convert(Vector{Float64}, peo[:,vind[3]]) convert(Vector{Float64}, peo[:,vind[4]]) convert(Vector{Float64}, peo[:,vind[5]]) convert(Vector{Float64}, peo[:,vind[6]])]*modelparameters' + rand(d, siz)
  mat7 = [ convert(Vector{Float64}, peo[:,viind[1]]) convert(Vector{Float64}, peo[:,viind[2]]) convert(Vector{Float64}, peo[:,viind[3]]) convert(Vector{Float64}, peo[:,viind[4]]) convert(Vector{Float64}, peo[:,viind[5]]) convert(Vector{Float64}, peo[:,viind[6]])]*modelparameters' + rand(d, siz)
  mat8 = [ convert(Vector{Float64}, peo[:,viiind[1]]) convert(Vector{Float64}, peo[:,viiind[2]]) convert(Vector{Float64}, peo[:,viiind[3]]) convert(Vector{Float64}, peo[:,viiind[4]]) convert(Vector{Float64}, peo[:,viiind[5]]) convert(Vector{Float64}, peo[:,viiind[6]])]*modelparameters' + rand(d, siz)
  mat9 = [ convert(Vector{Float64}, peo[:,ixnd[1]]) convert(Vector{Float64}, peo[:,ixnd[2]]) convert(Vector{Float64}, peo[:,ixnd[3]]) convert(Vector{Float64}, peo[:,ixnd[4]]) convert(Vector{Float64}, peo[:,ixnd[5]]) convert(Vector{Float64}, peo[:,ixnd[6]])]*modelparameters' + rand(d, siz)
  mat10 = [ convert(Vector{Float64}, peo[:,xnd[1]]) convert(Vector{Float64}, peo[:,xnd[2]]) convert(Vector{Float64}, peo[:,xnd[3]]) convert(Vector{Float64}, peo[:,xnd[4]]) convert(Vector{Float64}, peo[:,xnd[5]]) convert(Vector{Float64}, peo[:,xnd[6]])]*modelparameters' + rand(d, siz)
  mat11 = [ convert(Vector{Float64}, peo[:,xind[1]]) convert(Vector{Float64}, peo[:,xind[2]]) convert(Vector{Float64}, peo[:,xind[3]]) convert(Vector{Float64}, peo[:,xind[4]]) convert(Vector{Float64}, peo[:,xind[5]]) convert(Vector{Float64}, peo[:,xind[6]])]*modelparameters' + rand(d, siz)
  if size(entrants, 2) > 1
    allfids = [convert(Vector{Int64}, peo[:,fidnd[1]]) convert(Vector{Int64}, peo[:,fidnd[2]]) convert(Vector{Int64}, peo[:,fidnd[3]]) convert(Vector{Int64}, peo[:,fidnd[4]]) convert(Vector{Int64}, peo[:,fidnd[5]]) convert(Vector{Int64}, peo[:,fidnd[6]]) convert(Vector{Int64}, peo[:,fidnd[7]]) convert(Vector{Int64}, peo[:,fidnd[8]]) convert(Vector{Int64}, peo[:,fidnd[9]]) convert(Vector{Int64}, peo[:,fidnd[10]]) convert(Vector{Int64}, peo[:,fidnd[11]]) repmat(entfids, siz, 1)]
    entutil = EntrantsU(peo, entrants, modelparameters)
    for i = 1:size(peo, 1)
    #  outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i] entutil[i]])]
       best = indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i] entutil[i,1]])
       if best <= 11
         outp[i] = allfids[i, best]
       else
         outp[i] = allfids[i,11 + convert(Int,entutil[i,2])]
       end
    end
  else # there are no entrants, so "entrants" above is [0.0]', which has size(entrants, 2) == 1
    allfids = [convert(Vector{Int64}, peo[:,fidnd[1]]) convert(Vector{Int64}, peo[:,fidnd[2]]) convert(Vector{Int64}, peo[:,fidnd[3]]) convert(Vector{Int64}, peo[:,fidnd[4]]) convert(Vector{Int64}, peo[:,fidnd[5]]) convert(Vector{Int64}, peo[:,fidnd[6]]) convert(Vector{Int64}, peo[:,fidnd[7]]) convert(Vector{Int64}, peo[:,fidnd[8]]) convert(Vector{Int64}, peo[:,fidnd[9]]) convert(Vector{Int64}, peo[:,fidnd[10]]) convert(Vector{Int64}, peo[:,fidnd[11]])]
    for i = 1:size(peo, 1)
      outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i]])]
    end
  end
  return outp
end


# use countmap(choicemade) to count the results (!)  So easy.



#allents = zeros(size(people)[1], convert(Int, maxfid + floor((size(entrants)[2])/6)))
#=
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

=#


#=

# This version does not pass the distribution.  It is no slower than the verion which *does* pass the distribution
# 0.576737 seconds (8.93 M allocations: 329.331 MB)
function Demand(peo::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), persloc = [183 184] , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162] )
# constants/outputs/setup
  outp = zeros(siz)
  entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
# Computed utilities + error
  mat1 = [ convert(Vector{Float64}, peo[:,ind[1]]) convert(Vector{Float64}, peo[:,ind[2]]) convert(Vector{Float64}, peo[:,ind[3]]) convert(Vector{Float64}, peo[:,ind[4]]) convert(Vector{Float64}, peo[:,ind[5]]) convert(Vector{Float64}, peo[:,ind[6]])]*modelparameters' + rand(d, siz)
  mat2 = [ convert(Vector{Float64}, peo[:,iind[1]]) convert(Vector{Float64}, peo[:,iind[2]]) convert(Vector{Float64}, peo[:,iind[3]]) convert(Vector{Float64}, peo[:,iind[4]]) convert(Vector{Float64}, peo[:,iind[5]]) convert(Vector{Float64}, peo[:,iind[6]])]*modelparameters' + rand(d, siz)
  mat3 = [ convert(Vector{Float64}, peo[:,iiind[1]]) convert(Vector{Float64}, peo[:,iiind[2]]) convert(Vector{Float64}, peo[:,iiind[3]]) convert(Vector{Float64}, peo[:,iiind[4]]) convert(Vector{Float64}, peo[:,iiind[5]]) convert(Vector{Float64}, peo[:,iiind[6]])]*modelparameters' + rand(d, siz)
  mat4 = [ convert(Vector{Float64}, peo[:,ivnd[1]]) convert(Vector{Float64}, peo[:,ivnd[2]]) convert(Vector{Float64}, peo[:,ivnd[3]]) convert(Vector{Float64}, peo[:,ivnd[4]]) convert(Vector{Float64}, peo[:,ivnd[5]]) convert(Vector{Float64}, peo[:,ivnd[6]])]*modelparameters' + rand(d, siz)
  mat5 = [ convert(Vector{Float64}, peo[:,vnd[1]]) convert(Vector{Float64}, peo[:,vnd[2]]) convert(Vector{Float64}, peo[:,vnd[3]]) convert(Vector{Float64}, peo[:,vnd[4]]) convert(Vector{Float64}, peo[:,vnd[5]]) convert(Vector{Float64}, peo[:,vnd[6]])]*modelparameters' + rand(d, siz)
  mat6 = [ convert(Vector{Float64}, peo[:,vind[1]]) convert(Vector{Float64}, peo[:,vind[2]]) convert(Vector{Float64}, peo[:,vind[3]]) convert(Vector{Float64}, peo[:,vind[4]]) convert(Vector{Float64}, peo[:,vind[5]]) convert(Vector{Float64}, peo[:,vind[6]])]*modelparameters' + rand(d, siz)
  mat7 = [ convert(Vector{Float64}, peo[:,viind[1]]) convert(Vector{Float64}, peo[:,viind[2]]) convert(Vector{Float64}, peo[:,viind[3]]) convert(Vector{Float64}, peo[:,viind[4]]) convert(Vector{Float64}, peo[:,viind[5]]) convert(Vector{Float64}, peo[:,viind[6]])]*modelparameters' + rand(d, siz)
  mat8 = [ convert(Vector{Float64}, peo[:,viiind[1]]) convert(Vector{Float64}, peo[:,viiind[2]]) convert(Vector{Float64}, peo[:,viiind[3]]) convert(Vector{Float64}, peo[:,viiind[4]]) convert(Vector{Float64}, peo[:,viiind[5]]) convert(Vector{Float64}, peo[:,viiind[6]])]*modelparameters' + rand(d, siz)
  mat9 = [ convert(Vector{Float64}, peo[:,ixnd[1]]) convert(Vector{Float64}, peo[:,ixnd[2]]) convert(Vector{Float64}, peo[:,ixnd[3]]) convert(Vector{Float64}, peo[:,ixnd[4]]) convert(Vector{Float64}, peo[:,ixnd[5]]) convert(Vector{Float64}, peo[:,ixnd[6]])]*modelparameters' + rand(d, siz)
  mat10 = [ convert(Vector{Float64}, peo[:,xnd[1]]) convert(Vector{Float64}, peo[:,xnd[2]]) convert(Vector{Float64}, peo[:,xnd[3]]) convert(Vector{Float64}, peo[:,xnd[4]]) convert(Vector{Float64}, peo[:,xnd[5]]) convert(Vector{Float64}, peo[:,xnd[6]])]*modelparameters' + rand(d, siz)
  mat11 = [ convert(Vector{Float64}, peo[:,xind[1]]) convert(Vector{Float64}, peo[:,xind[2]]) convert(Vector{Float64}, peo[:,xind[3]]) convert(Vector{Float64}, peo[:,xind[4]]) convert(Vector{Float64}, peo[:,xind[5]]) convert(Vector{Float64}, peo[:,xind[6]])]*modelparameters' + rand(d, siz)
  if size(entrants, 2) > 1
    allfids = [convert(Vector{Int64}, peo[:,fidnd[1]]) convert(Vector{Int64}, peo[:,fidnd[2]]) convert(Vector{Int64}, peo[:,fidnd[3]]) convert(Vector{Int64}, peo[:,fidnd[4]]) convert(Vector{Int64}, peo[:,fidnd[5]]) convert(Vector{Int64}, peo[:,fidnd[6]]) convert(Vector{Int64}, peo[:,fidnd[7]]) convert(Vector{Int64}, peo[:,fidnd[8]]) convert(Vector{Int64}, peo[:,fidnd[9]]) convert(Vector{Int64}, peo[:,fidnd[10]]) convert(Vector{Int64}, peo[:,fidnd[11]]) repmat(entfids, siz, 1)]
    entutil = EntrantsU(peo, entrants, modelparameters)
    for i = 1:size(peo, 1)
    #  outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i] entutil[i]])]
       best = indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i] entutil[i,1]])
       if best <= 11
         outp[i] = allfids[i, best]
       else
         outp[i] = allfids[i,11 + convert(Int,entutil[i,2])]
       end
    end
  else # there are no entrants, so "entrants" above is [0.0]', which has size(entrants, 2) == 1
    allfids = [convert(Vector{Int64}, peo[:,fidnd[1]]) convert(Vector{Int64}, peo[:,fidnd[2]]) convert(Vector{Int64}, peo[:,fidnd[3]]) convert(Vector{Int64}, peo[:,fidnd[4]]) convert(Vector{Int64}, peo[:,fidnd[5]]) convert(Vector{Int64}, peo[:,fidnd[6]]) convert(Vector{Int64}, peo[:,fidnd[7]]) convert(Vector{Int64}, peo[:,fidnd[8]]) convert(Vector{Int64}, peo[:,fidnd[9]]) convert(Vector{Int64}, peo[:,fidnd[10]]) convert(Vector{Int64}, peo[:,fidnd[11]])]
    for i = 1:size(peo, 1)
      outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i]])]
    end
  end
  return outp
end


=#

















#end
