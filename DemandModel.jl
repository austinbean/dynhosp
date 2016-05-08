# Use the parameters of the demand model to estimate choices:

# Version: 05 07 16

using DataFrames
using DataArrays
using Distributions

# NEED TO ADD BEDS HERE -
people = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);

# Check the NFP status variable in the above

modcoeffs = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Model.csv", header = true);
distance_c = modcoeffs[1, 2]
distsq_c = modcoeffs[2, 2]
neoint_c = modcoeffs[3, 2]
soloint_c = modcoeffs[4, 2]
closest_c = modcoeffs[5, 2]
distbed_c = modcoeffs[6, 2]
# ADD BEDS HERE VIA STATA


# This is needed to clean out the missing values among fids.  Changes them to 0.
maxfid = 11;
for i in names(people)
  # This won't work for string variables, obviously.  Fix later. 
  people[isna(people[i]), i] = 0
  # eval(ex1)
end



function fidfinder(fidvect::Array{Int64, 2}, choices::DataFrame; maxfid = 11)
    #=
      This function generates a vector of Booleans which index the rows in the
      dataframe in which the individual has some hospital with fid in fidvect
      as an option.  It starts with all falses and iteratively takes subset | (result)
      which will be true when the result expression is true.  It operates on a whole DataFrame
    =#
    subset = falses(size(choices)[1])
      for k in 1:size(fidvect)[1]
        targ = fidvect[k]
        for j = 1:maxfid
           subex = parse("people[:fid$j].==$targ")
           true_v = eval(subex)
      #     print( size(true_v), "   ") # for testing purposes
           subset = subset | true_v
        end
      end
    return subset
end

# The next function will find values in the one-row dataframe element given a list of symbols

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
    return RowSizeError
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
        el = convert(Int64, el)
        (loc, symb) = rowfindfid(choicerow, el)
        if loc != 0
          fid_num = replace(string(symb), r"[fid]", "")
          # Change NeoIntensive in the choice row to the value in the state row
          neo = Symbol("NeoIntensive"*fid_num)
          choicerow[neo] = convert(Int64, staterow[findfirst(staterow, el)+2])
          # Change SoloIntermeidate in the choice row to the value in the state row
          solo = Symbol("SoloIntermediate"*fid_num)
          choicerow[solo] = convert(Int64, staterow[findfirst(staterow, el)+1])
        else
          print(el, " not found in row ")
        end
      end
    end
end


dist_μ = 0;
dist_σ = 1;
dist_ξ = 0;
srand(123)
d = GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# I do need the constant:
γ = eulergamma;

# We need to do several things:
#=
 - Remove all NA values from the people frame and replace them with 0's    ✓
 - Write some function (variable argument numbers?) which finds a vector of fids in all people spots ✓

=#


modelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]

function DemandModel(individuals::DataFrame, modelparameters::Array{Float64, 2})
  #=
    The goal for this function is to -
      - take the whole set of people, compute the deterministic components of utility, add the random shock, find the maximizer
      - count the number maximized by fid: this will be the demand.
  =#



end
