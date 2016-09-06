# Use the parameters of the demand model to estimate choices:

# Version: 05 09 16
#
# using DataFrames
# using DataArrays
# using Distributions
#

# Now in use in Main.jl

# Indices:
#people = DataFrames.readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 1 Individual Choices.csv", header = true);
# people = DataFrames.readtable(pathpeople*"TX 2005 Individual Choices.csv", header = true);
#
# fid1loc = people.colindex.lookup[:fid1]
# fid2loc = people.colindex.lookup[:fid2]
# fid3loc = people.colindex.lookup[:fid3]
# fid4loc = people.colindex.lookup[:fid4]
# fid5loc = people.colindex.lookup[:fid5]
# fid6loc = people.colindex.lookup[:fid6]
# fid7loc = people.colindex.lookup[:fid7]
# fid8loc = people.colindex.lookup[:fid8]
# fid9loc = people.colindex.lookup[:fid9]
# fid10loc = people.colindex.lookup[:fid10]
# fid11loc = people.colindex.lookup[:fid11]
# fid1loc = 2; fid2loc = 11; fid3loc = 20; fid4loc = 29; fid5loc = 38; fid6loc = 47; fid7loc = 56; fid8loc = 65; fid9loc = 74; fid10loc = 83; fid11loc = 92
# fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92]
#fid1loc = 2, fid2loc = 11, fid3loc = 20, fid4loc = 29, fid5loc = 38, fid6loc = 47, fid7loc = 56, fid8loc = 65, fid9loc = 74, fid10loc = 83, fid11loc = 92

function fidfinder(fidvect::Array{Int64, 2}, choices::Matrix; maxfid = 11, fid1loc = 2, fid2loc = 11, fid3loc = 20, fid4loc = 29, fid5loc = 38, fid6loc = 47, fid7loc = 56, fid8loc = 65, fid9loc = 74, fid10loc = 83, fid11loc = 92)
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

function fidfinder(fidvect::Array{Int64, 1}, choices::Matrix; maxfid = 11, fid1loc = 2, fid2loc = 11, fid3loc = 20, fid4loc = 29, fid5loc = 38, fid6loc = 47, fid7loc = 56, fid8loc = 65, fid9loc = 74, fid10loc = 83, fid11loc = 92)
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





# fid1loc = 2, fid2loc = 18, fid3loc = 34, fid4loc = 50, fid5loc = 66, fid6loc = 82, fid7loc = 98, fid8loc = 114, fid9loc = 130, fid10loc = 146, fid11loc = 162
# fidnd = [2, 18, 34, 50, 66, 82, 98, 114, 130, 146, 162]

# hisrow1 = [1131021 -1 -1 0 1 1 1 9 9 9 9.0]
# mfids = [1131021]
# state history/hisrow  has the form: [ fid, act_solo, act_int, choice prob, action taken, demand realized, perturbed] × (# facilities)  ⋃ [ level1's total, level2's total, level3's total, aggregate prob]
# choicerow/people  has the form: [identity, fid, facility, Total Beds, NeoIntensive, TotalDeliveries, Transfers Out No NICU, Transfers In Has NICU, Transfers Out Has NICU, Not For Profit Status (#), Solo Intermediate, distance, Is Closest?, Selected?, NFP ?, distance × bed, distance²] × (# facilities) ⋃ [Patient Zip, CMS MDC, APR MDC, CMS DRG, APR DRG, Zip Lat, Zip Long]
# new: [fid, NeoIntensive, Solo Intermediate, distance, Is Closest?, Selected?, distance × bed, distance², amount charged]
# @time: 0.070399 seconds (1.60 M allocations: 66.686 MB, 15.01% gc time) - this for one value in mfids.
# @time: 0.148658 seconds (2.17 M allocations: 88.855 MB) - for 10 values in mfids.  choiceintloc was 3, choicesololoc was 9
function rowchange(hisrow::Array{Float64, 2}, mfids::Array{Int64}, people::Matrix; choiceintloc = 1, choicesololoc = 2, lenrow = (maximum(size(hisrow))-4), fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92]) # , hisfd = collect(1:7:lenrow)
  for i in 1:size(people, 1)
    change_fids = intersect(slice(people, i, fidnd), mfids) # 14 allocations: 464 bytes
    for j in change_fids
      frm  =  findfirst(hisrow ,j)  # 5 allocations / 208 bytes
      to =  findfirst(slice(people, i, :),j)  # 12 allocations / 352 bytes
      people[i, to + choicesololoc] = hisrow[frm + 1] # 5 allocations / 176 bytes
      people[i, to + choiceintloc] =  hisrow[frm + 2] # 5 allocations / 176 bytes
    end
  end # when someone exits, can you not search for that fid?
  return people
end

# sample entrants1 = [99999 1 0 120 32.96  -96.8385] [newrow[fidloc] newrow[act_intloc] newrow[act_sololoc] entrantbeds ent_lat ent_lon]
# sample entrants2 = [99999 1 0 120 32.96  -96.8385 888888 0 1 120 32.96  -96.8385]

function EntrantsU(peo::Matrix, entrants::Array{Float64, 2}, modelparameters::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ), persloc = [ 105 106], entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize))
  siz = size(peo,1)
  rands = rand(d, siz, entnum)
  entvals = zeros(siz, entnum)
  plocs = peo[:,persloc[:]] # this format is: LATITUDE, LONGITUDE
  for j = 1:entnum
    for i = 1:siz # ENTRANT FORMAT ends with [Latitude, Longitude]
      d1 = distance(entrants[6*j-1], entrants[6*j], plocs[i,1], plocs[i,2])
      if d1 < 50
        entvals[i,j] += d1*modelparameters[1]
        entvals[i,j] += ((d1)^2)*modelparameters[2]
        entvals[i,j] += entrants[6*j-3]*modelparameters[3]
        entvals[i,j] += entrants[6*j-4]*modelparameters[4]
    #    entvals[i,j] += 0*modelparameters[5] # this is specifying that "closest" is always 0 for entrants.  It can be fixed, but would be really annoying.
        entvals[i,j] += ((d1)*entrants[6*j-2]/100)*modelparameters[6]
      else
        entvals[i,j] = -99 # set value to large negative number when distance is too large: won't be chosen
      end
    end
  end
  utils, inds = findmax(entvals + rands, 2) # findmax returns value and index over given dimension
  return [utils ind2sub(size(entvals), vec(inds))[2]]  # note that due to the randomization, this will generally not return -999, but -999 + rand
end

# Another version with fewer allocations:
# entrants0 = Array{Float64, 2}()
# entrants1 = [99999 1 0 120 32.96  -96.8385]
# Much improved now.
# Most of the slowness comes when entrants are present.
# choicerow/people  has the form: [identity] ∪ [fid, NeoIntensive, Solo Intermediate, distance, Is Closest?, Selected?, distance × bed, distance², amount charged] × (# facilities) ⋃ [Patient Zip, DRG, medicaid, Private insurance, Zip Lat, Zip Long]
# Order of Demand Model Coefficients: [distance_c  distsq_c  neoint_c  soloint_c  closest_c  distbed_c ]
# demandmodelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]
# Original locations, pre breaking change: ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162]
# Prior to changing return type to include zip.
# @time DemandModel(pinsured, privatedemandmodelparameters, entrants)
#  0.527198 seconds (2.45 M allocations: 272.934 MB, 25.48% gc time) - but this was for a larger set, probably mothers and infants.
#ind = [5 9 3 4 6 8 ]; iind = [14 18 12 13 15 17 ]; iiind = [23 27 21 22 24 26 ]; ivnd = [32 36 30 31 33 35 ]; vnd = [41 45 39 40 42 44 ]; vind = [50 54 48 49 51 53 ]; viind = [59 63 57 58 60 62 ]; viiind = [68 72 66 67 69 71 ]; ixnd = [77 81 75 76 78 80 ]; xnd = [86 90 84 85 87 89 ]; xind = [95 99 93 94 96 98 ]

function DemandModel(peo::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ), ziploc = 101, drgloc = 102, entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92] , ind = [5 9 3 4 6 8 ], iind = [14 18 12 13 15 17 ], iiind = [23 27 21 22 24 26 ], ivnd = [32 36 30 31 33 35 ], vnd = [41 45 39 40 42 44 ], vind = [50 54 48 49 51 53 ], viind = [59 63 57 58 60 62 ], viiind = [68 72 66 67 69 71 ], ixnd = [77 81 75 76 78 80 ], xnd = [86 90 84 85 87 89 ], xind = [95 99 93 94 96 98 ] )
# Computed utilities + error
  rand_el = Array(Float64, siz)
  mat1 = peo[:,ind[1:6]]*modelparameters' + rand!(d, rand_el)
  mat2 = peo[:,iind[1:6]]*modelparameters' + rand!(d, rand_el)
  mat3 = peo[:,iiind[1:6]]*modelparameters' + rand!(d, rand_el)
  mat4 = peo[:,ivnd[1:6]]*modelparameters' + rand!(d, rand_el)
  mat5 = peo[:,vnd[1:6]]*modelparameters' + rand!(d, rand_el)
  mat6 = peo[:,vind[1:6]]*modelparameters' + rand!(d, rand_el)
  mat7 = peo[:,viind[1:6]]*modelparameters' + rand!(d, rand_el)
  mat8 = peo[:,viiind[1:6]]*modelparameters' + rand!(d, rand_el)
  mat9 = peo[:,ixnd[1:6]]*modelparameters' + rand!(d, rand_el)
  mat10 = peo[:,xnd[1:6]]*modelparameters' + rand!(d, rand_el)
  mat11 = peo[:,xind[1:6]]*modelparameters' + rand!(d, rand_el)
  if size(entrants, 2) > 1
    entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
    allfids = [peo[:,fidnd[1:11]] repmat(entfids, siz, 1)] #maybe fill(entfids, siz) would work faster?  Maybe speed is the same but allocations lower.
    entutil = EntrantsU(peo, entrants, modelparameters)
    vals, inds = findmax([mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11 entutil[:,1]] , 2)
    # Fuck-up is right here.
    outp = map( (i,x)->allfids[i,x], collect(1:size(mat1,1)), ind2sub((size(mat1,1),12), vec(inds) )[2] )
  else #  no entrants
      allfids = peo[:, fidnd[1:11]]
      vals, inds = findmax([mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11] , 2)
      outp = map((i,x)->allfids[i,x], 1:size(mat1,1), ind2sub((size(mat1,1),11), vec(inds) )[2] )
  end
return hcat( peo[:, ziploc], peo[:, drgloc], outp)
end

# Call DetUtil first, then this.
function DemandModel2(detutil::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ), ziploc = 1, drgloc = 2, entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(detutil,1), fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92], ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
# Computed utilities + error
  rand_el = Array{Float64}(siz, 11)
  if size(entrants, 2) > 1
    entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
    entutil = EntrantsU(XXXdetutil??XXX , entrants, modelparameters)
    vals, inds = findmax( [detutil[:,ulocs[:]] + rand!(d, rand_el) entutil[:,1]] , 2)
    # Fuck-up is right here.
    outp = map( (i,x)->allfids[i,x], collect(1:size(mat1,1)), ind2sub((size(mat1,1),12), vec(inds) )[2] )
  else #  no entrants
      vals, inds = findmax(detutil[:,ulocs[:]] + rand!(d, rand_el), 2) # returns indices in the range [1, ..., 11]
      outp = map((i,x)->detutil[i,x], 1:siz, 2*(ind2sub((siz,11), vec(inds) )[2])+2 )
  end
return hcat( detutil[:, ziploc], detutil[:, drgloc], outp)
end



#=
# This next loop will list the column index of the relevant parameters in order.
for i = 1:11
  println("Group ", i)
  print("[")
  for name = ["distance"  "distsq" "NeoIntensive" "SoloIntermediate" "closest" "dist_bed" "amountcharged"] # can add "amountcharged"
    print(people.colindex.lookup[convert(Symbol, name*string(i))], " " )
  end
  print("]")
  println("    ")
end
=#

# entrants = Array{Float64,2}(); dist_μ = 0; dist_σ = 1; dist_ξ = 0; d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ); entsize = 6; entnum = convert(Int, size(entrants, 2)/entsize); siz = size(peo,1); persloc = [183 184]; ind = [12 17 11 5 13 16]; iind = [28 33 27 21 29 32]; iiind = [44 49 43 37 45 48]; ivnd = [60 65 59 53 61 64]; vnd = [76 81 75 69 77 80]; vind = [92 97 91 85 93 96]; viind = [108 113 107 101 109 112]; viiind = [124 129 123 117 125 128]; ixnd = [140 145 139 133 141 144]; xnd = [156 161 155 149 157 160]; xind = [172 177 171 165 173 176]; fidnd = [2 18 34 50 66 82 98 114 130 146 162]
#=
function DM2(peo::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ), entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), persloc = [183 184] , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162] )
# Computed utilities + error
# change_fids = intersect(slice(people, i, fidnd), mfids)
  rand_el = Array(Float64, siz)
  vals = fill(typemin(Float64), siz)
  inds = fill(0, siz)
  for indices in [ind, iind, iiind, ivnd, vnd, vind, viind, viiind, ixnd, xnd, xind]
    mat = peo[:, indices[1:6]]*modelparameters'+rand!(d, rand_el)

  end
  if size(entrants, 2) > 1
    entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
    allfids = [peo[:,fidnd[1:11]] repmat(entfids, siz, 1)]
    entutil = EntrantsU(peo, entrants, modelparameters)
    vals, inds = findmax([mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11 entutil[:,1]] , 2)
    outp = map( (i,x)->allfids[i,x], collect(1:size(mat1,1)), ind2sub(( size(mat1,1),12 ), vec(inds) )[2] )
    outp = [ allfids[i,x] for i = 1:siz, x = ind2sub((siz,12), vec(inds))[2] ]
  else #  no entrants
      allfids = peo[:, fidnd[1:11]]
      vals, inds = findmax([mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11] , 2)
      outp = map((i,x)->allfids[i,x], 1:size(mat1,1), ind2sub((size(mat1,1),11), vec(inds) )[2] )
  end
return outp
end
=#


##### Tester:
#=
function DemandMod3(peo::Array{Float64,2}, modelparameters::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ),  siz = size(peo,1), persloc = [183 184] , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162] )
  mat1 = peo[:,ind[1:6]]*modelparameters' + rand(d, siz)
  mat2 = peo[:,iind[1:6]]*modelparameters' + rand(d, siz)
  mat3 = peo[:,iiind[1:6]]*modelparameters' + rand(d, siz)
  mat4 = peo[:,ivnd[1:6]]*modelparameters' + rand(d, siz)
  mat5 = peo[:,vnd[1:6]]*modelparameters' + rand(d, siz)
  mat6 = peo[:,vind[1:6]]*modelparameters' + rand(d, siz)
  mat7 = peo[:,viind[1:6]]*modelparameters' + rand(d, siz)
  mat8 = peo[:,viiind[1:6]]*modelparameters' + rand(d, siz)
  mat9 = peo[:,ixnd[1:6]]*modelparameters' + rand(d, siz)
  mat10 = peo[:,xnd[1:6]]*modelparameters' + rand(d, siz)
  mat11 = peo[:,xind[1:6]]*modelparameters' + rand(d, siz)
  allfids = peo[:, fidnd[1:11]]
  vals, inds = findmax([mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11] , 2)
  outp = map((i,x)->allfids[i,x], collect(1:size(mat1,1)), ind2sub((size(mat1,1),12), vec(inds) )[2] )
  return outp
end
=#
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





#=

# 1.599738 seconds (9.47 M allocations: 355.385 MB, 56.89% gc time)
# 0.533586 seconds (7.71 M allocations: 310.720 MB)
# w/ 400,000 rows:  0.518312 seconds (4.55 M allocations: 601.591 MB) (removing the creation of allfids at the bottom when there are no entrants)
# entrants = Array{Float64, 2}() # for a fast test with no entrants
# That's more costly than I would like
# Experiment with eliminating the conversions, just form the matrices.
function DemandModel(peo::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; dist_μ = 0, dist_σ = 1, dist_ξ = 0, d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ), entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), persloc = [183 184] , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162] )
#  d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)
# constants/outputs/setup
  outp = zeros(siz)
  entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
# Computed utilities + error
  mat1 = [ peo[:,ind[1]] peo[:,ind[2]] peo[:,ind[3]]  peo[:,ind[4]] peo[:,ind[5]] peo[:,ind[6]]]*modelparameters' + rand(d, siz)
  mat2 = [ peo[:,iind[1]] peo[:,iind[2]] peo[:,iind[3]]  peo[:,iind[4]] peo[:,iind[5]] peo[:,iind[6]]]*modelparameters' + rand(d, siz)
  mat3 = [ peo[:,iiind[1]] peo[:,iiind[2]] peo[:,iiind[3]]  peo[:,iiind[4]] peo[:,iiind[5]] peo[:,iiind[6]]]*modelparameters' + rand(d, siz)
  mat4 = [ peo[:,ivnd[1]] peo[:,ivnd[2]] peo[:,ivnd[3]]  peo[:,ivnd[4]] peo[:,ivnd[5]] peo[:,ivnd[6]]]*modelparameters' + rand(d, siz)
  mat5 = [ peo[:,vnd[1]] peo[:,vnd[2]] peo[:,vnd[3]]  peo[:,vnd[4]] peo[:,vnd[5]] peo[:,vnd[6]]]*modelparameters' + rand(d, siz)
  mat6 = [ peo[:,vind[1]] peo[:,vind[2]] peo[:,vind[3]]  peo[:,vind[4]] peo[:,vind[5]] peo[:,vind[6]]]*modelparameters' + rand(d, siz)
  mat7 = [ peo[:,viind[1]] peo[:,viind[2]] peo[:,viind[3]]  peo[:,viind[4]] peo[:,viind[5]] peo[:,viind[6]]]*modelparameters' + rand(d, siz)
  mat8 = [ peo[:,viiind[1]] peo[:,viiind[2]] peo[:,viiind[3]]  peo[:,viiind[4]] peo[:,viiind[5]] peo[:,viiind[6]]]*modelparameters' + rand(d, siz)
  mat9 = [ peo[:,ixnd[1]] peo[:,ixnd[2]] peo[:,ixnd[3]]  peo[:,ixnd[4]] peo[:,ixnd[5]] peo[:,ixnd[6]]]*modelparameters' + rand(d, siz)
  mat10 = [ peo[:,xnd[1]] peo[:,xnd[2]] peo[:,xnd[3]]  peo[:,xnd[4]] peo[:,xnd[5]] peo[:,xnd[6]]]*modelparameters' + rand(d, siz)
  mat11 = [ peo[:,xind[1]] peo[:,xind[2]] peo[:,xind[3]]  peo[:,xind[4]] peo[:,xind[5]] peo[:,xind[6]]]*modelparameters' + rand(d, siz)
  if size(entrants, 2) > 1
    allfids = [peo[:,fidnd[1]] peo[:,fidnd[2]] peo[:,fidnd[3]] peo[:,fidnd[4]] peo[:,fidnd[5]] peo[:,fidnd[6]] peo[:,fidnd[7]] peo[:,fidnd[8]] peo[:,fidnd[9]] peo[:,fidnd[10]] peo[:,fidnd[11]] repmat(entfids, siz, 1)]
    entutil = EntrantsU(peo, entrants, modelparameters)
    for i = 1:size(peo, 1)
    #  outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i] entutil[i]])]
       best = indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i] entutil[i,1]])
       # mat = [mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11 entutil[:,1]]
       # vals, inds = findmax([mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11 entutil[:,1]] , 2)
       # ind2sub(size(mat), vec(inds) )[2]
       # map( (i,x)->allfids[i,x], collect(1:size(mat,1)), ind2sub(size(mat), vec(inds) )[2] )
       if best <= 11
         outp[i] = allfids[i, best]
       else
         outp[i] = allfids[i,11 + convert(Int,entutil[i,2])]
       end
    end
  else # there are no entrants, so "entrants" above is [0.0]', which has size(entrants, 2) == 1
  #  allfids = [ peo[:,fidnd[1]] peo[:,fidnd[2]] peo[:,fidnd[3]] peo[:,fidnd[4]] peo[:,fidnd[5]] peo[:,fidnd[6]] peo[:,fidnd[7]] peo[:,fidnd[8]] peo[:,fidnd[9]] peo[:,fidnd[10]] peo[:,fidnd[11]]]
    for i = 1:size(peo, 1)
      # outp[i] = allfids[i, indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i]])]
      outp[i] = peo[i, fidnd[indmax([mat1[i] mat2[i] mat3[i] mat4[i] mat5[i] mat6[i] mat7[i] mat8[i] mat9[i] mat10[i] mat11[i]])]]
    end
  end
  return outp
end
=#



# DEPRECATED:
# function rowfindfid(targ::DataFrame, value::Int64; vals = [:fid1, :fid2, :fid3, :fid4, :fid5, :fid6, :fid7, :fid8, :fid9, :fid10, :fid11] )
#   #=
#     This function takes a single row of a dataframe "targ" and finds the value in "value"
#     by looking in the columns in "vals".  Those items in vals must be a list of symbols, so
#     [ :fid1, :fid2, ..., :fidN ] (commas necessary).  If the DataFrame is too big, it Returns
#     the RowSizeError exception defined above.  To search specific fids, must call function as
#     rowfindfid(targ, value, vals = [:fidx, :fidy])
#   =#
#   if size(targ)[1] > 1
#     return RowSizeError #defined in Main.jl
#   end
#   numb = 0
#   index = :ident
#   for i in vals # what about NAs?  All of the fid columns are zeros.
#       if targ[i][1] == value # the second [1] is needed to access the *value* of the element, rather than the 1 element DataFrame
#           numb = targ.colindex.lookup[i] # returns the *column number* of the element i
#           index =  i
#       end
#   end
#   return numb, index
# end







#end
