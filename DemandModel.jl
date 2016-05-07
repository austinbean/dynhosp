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
    subset = falses(size(choices)[1])
      for k in 1:size(fidvect)[1]
        targ = fidvect[k]
        for j = 1:maxfid
           subex = parse("choices[:fid$j].==$targ")
           true_v = eval(subex)
      #     print( size(true_v), "   ") # for testing purposes
           subset = subset | true_v
        end
      end  
    return subset
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
