# Use the parameters of the demand model to estimate choices:

using DataFrames
using DataArrays
using Distributions


people = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);

# Check the NFP status variable in the above

modcoeffs = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Model.csv", header = true);


# This is needed to clean out the missing values among fids.  Changes them to 0.
maxfid = 11;
for i = 1:maxfid
  ex1 = parse("people[isna(people[:fid$i]), :fid$i] = 0")
  eval(ex1)
end



function fidfinder(fidvect::Array{Int64, 2}, choices::DataFrame; maxfid = 11)
    #=
      This function generates a vector of Booleans which index the rows in the
      dataframe in which the individual has some hospital with fid in fidvect
      as an option.
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


function rowchange(staterow::Array{Float64,2}, choicerow::DataFrame; endfields_state = 4, fields_state = 7, fields_people = 15, endfields_people = 7)
  #=
     This function should take a row of the state history (staterow), and a row of
     the choices (choicerow) and:
     1.  determines the number of fids in the staterow
     2.  Determines the number of fids in the choicerow
     3.  When a fid in the staterow matches a fid in the choicerow, map the values
         from the staterow to the choicerow
   Notes - need to do something special for entrants.  Can just check if fid sets are overlapping, I guess.
   Once I know this, I also need to check whether the new hospital is the closest.

   staterow has the form: [ fid, act_solo, act_int, choice prob, action taken, demand realized, perturbed] × (# facilities)  ⋃ [ level1's total, level2's total, level3's total, aggregate prob]
   choicerow has the form: [identity, fid, facility, NeoIntensive, TotalDeliveries, Transfers Out No NICU, Transfers In Has NICU, Not For Profit Status (#), Solo Intermediate, distance, Is Closest?, Selected?, NFP ?, distance × bed, distance²] × (# facilities) ⋃ [Patient Zip, CMS MDC, APR MDC, CMS DRG, APR DRG, Zip Lat, Zip Long]

  =#
    mktnumfids = convert(Int, ((size(staterow)[2])-endfields_state)/fields_state) # number of facilities
    mktfids = [ el for el in staterow[1,1:fields_state:end-endfields_state]] # Collects the fids in the market

    peoplefids = [x for x in choicerow[1, 2:fields_people:end-endfields_people] ]

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
 - Write some function (variable argument numbers?) which finds a vector of fids in all people spots

=#


function DemandModel(individuals::DataFrame, modelparameters::Array{Float64, 2}, hospitalparameters::Array{Float64, 2})




end
