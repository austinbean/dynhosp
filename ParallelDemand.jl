# Trying to write demand as a function to be parallelized

# Try to get the demand computation to be as low on allocations and memory as possible.

using DataFrames
using DataArrays
using Distributions

include("/Users/austinbean/Desktop/dynhosp/Distance.jl")

people = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 Individual Choices.csv", header = true);

# Ok - what do I want?  Convert this to a matrix, then see if the multiplication and addition of
# the shock is faster and involves fewer allocations.

a = Set()
for el in people.columns
  push!(a, typeof(el))
end

for i in names(people)
  if ( typeof(people[i]) == DataArrays.DataArray{Float64,1} )
    people[isna(people[i]), i] = 0
  elseif (typeof(people[i]) == DataArrays.DataArray{Int64,1})
    people[isna(people[i]), i] = 0
  elseif typeof(people[i]) == DataArrays.DataArray{ByteString,1}
    people[isna(people[i]), i] = "NONE"
  end
  if sum(size(people[isna(people[i]), i]))>0
    print(i, "\n")
  end
end

modcoeffs = readtable("/Users/austinbean/Google Drive/Texas Inpatient Discharge/TX 2005 Model.csv", header = true);
distance_c = modcoeffs[1, 2]
distsq_c = modcoeffs[2, 2]
neoint_c = modcoeffs[3, 2]
soloint_c = modcoeffs[4, 2]
closest_c = modcoeffs[5, 2]
distbed_c = modcoeffs[6, 2]

modelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]

entrants = [0.0]'

dist_μ = 0;
dist_σ = 1;
dist_ξ = 0;
srand(123)
d = GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ)

type NotRowError <: Exception end



#=
Roadmap -

Benchmark time for current version of DemandModel in DemandModel.jl:
91.463204 seconds (130.02 M allocations: 2.696 GB, 1.30% gc time)

Current time for imptest32:
0.629812 seconds (9.09 M allocations: 575.473 MB, 13.98% gc time)


=#



function PDemand(people::DataFrame, frname::ASCIIString, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; maxfid = 11, ent_length = 6 )
  #=
    The goal for this function is to -
      - Write this function to be parallelizable:
        - operate on a row
        - return the value for the row.
      The first two arguments are a dataframe (people) and the NAME of that dataframe (frname) as a string.  This is important at the end
      - Entrants - measure length, computes number of entrants.
      - Format of entrants is: [fid act_int act_solo entrantbeds ent_lat ent_lon] x (# entrants)
      - Remember: distbed is distance * beds/100
      - Would it be faster to compute the sum in a different way?  Maybe the allocations result from computing each element of the
      sum separately?  That might be a speed up.
  =#
    distance_c  = modelparameters[1]
    distsq_c  = modelparameters[2]
    neoint_c  = modelparameters[3]
    soloint_c  = modelparameters[4]
    closest_c  = modelparameters[5]
    distbed_c = modelparameters[6]
    choicemade = zeros(size(people, 1))
      shock = rand(d, maxfid)
      if people[:fid1].data[1] != 0
        val1 = people[:distance1]*distance_c + people[:distsq1]*distsq_c + people[:SoloIntermediate1]*soloint_c + people[:NeoIntensive1]*neoint_c + people[:closest1]*closest_c + people[:dist_bed1]*distbed_c + shock[1]
        dist1 = people[:distance1]
      else
        val1 = -10^10
        dist1 = 999
      end
      if people[:fid2].data[1] != 0
        val2 = people[:distance2]*distance_c + people[:distsq2]*distsq_c + people[:SoloIntermediate2]*soloint_c + people[:NeoIntensive2]*neoint_c + people[:closest2]*closest_c + people[:dist_bed2]*distbed_c + shock[2]
        dist2 = people[:distance2]
      else
        val2 = -10^10
        dist2 = 999
      end
      if people[:fid3].data[1] != 0
        val3 = people[:distance3]*distance_c + people[:distsq3]*distsq_c + people[:SoloIntermediate3]*soloint_c + people[:NeoIntensive3]*neoint_c + people[:closest3]*closest_c + people[:dist_bed3]*distbed_c + shock[3]
        dist3 = people[:distance3]
      else
        val3 = -10^10
        dist3 = 999
      end
      if people[:fid4].data[1] != 0
        val4 = people[:distance4]*distance_c + people[:distsq4]*distsq_c + people[:SoloIntermediate4]*soloint_c + people[:NeoIntensive4]*neoint_c + people[:closest4]*closest_c + people[:dist_bed4]*distbed_c + shock[4]
        dist4 = people[:distance4]
      else
        val4 = -10^10
        dist4 = 999
      end
      if people[:fid5].data[1] != 0
        val5 = people[:distance5]*distance_c + people[:distsq5]*distsq_c + people[:SoloIntermediate5]*soloint_c + people[:NeoIntensive5]*neoint_c + people[:closest5]*closest_c + people[:dist_bed5]*distbed_c + shock[5]
        dist5 = people[:distance5]
      else
        val5 = -10^10
        dist5 = 999
      end
      if people[:fid6].data[1] != 0
        val6 = people[:distance6]*distance_c + people[:distsq6]*distsq_c + people[:SoloIntermediate6]*soloint_c + people[:NeoIntensive6]*neoint_c + people[:closest6]*closest_c + people[:dist_bed6]*distbed_c + shock[6]
        dist6 = people[:distance6]
      else
        val6 = -10^10
        dist6 = 999
      end
      if people[:fid7].data[1] != 0
        val7 = people[:distance7]*distance_c + people[:distsq7]*distsq_c + people[:SoloIntermediate7]*soloint_c + people[:NeoIntensive7]*neoint_c + people[:closest7]*closest_c + people[:dist_bed7]*distbed_c + shock[7]
        dist7 = people[:distance7]
      else
        val7 = -10^10
        dist7 = 999
      end
      if people[:fid8].data[1] != 0
        val8 = people[:distance8]*distance_c + people[:distsq8]*distsq_c + people[:SoloIntermediate8]*soloint_c + people[:NeoIntensive8]*neoint_c + people[:closest8]*closest_c + people[:dist_bed8]*distbed_c + shock[8]
        dist8 = people[:distance8]
      else
        val8 = -10^10
        dist8 = 999
      end
      if people[:fid9].data[1] != 0
        val9 = people[:distance9]*distance_c + people[:distsq9]*distsq_c + people[:SoloIntermediate9]*soloint_c + people[:NeoIntensive9]*neoint_c + people[:closest9]*closest_c + people[:dist_bed9]*distbed_c + shock[9]
        dist9 = people[:distance9]
      else
        val9 = -10^10
        dist9 = 999
      end
      if people[:fid10].data[1] != 0
        val10 = people[:distance10]*distance_c + people[:distsq10]*distsq_c + people[:SoloIntermediate10]*soloint_c + people[:NeoIntensive10]*neoint_c + people[:closest10]*closest_c + people[:dist_bed10]*distbed_c + shock[10]
        dist10 = people[:distance10]
      else
        val10 = -10^10
        dist10 = 999
      end
      if people[:fid11].data[1] != 0
        val11 = people[:distance11]*distance_c + people[:distsq11]*distsq_c + people[:SoloIntermediate11]*soloint_c + people[:NeoIntensive11]*neoint_c + people[:closest11]*closest_c + people[:dist_bed11]*distbed_c + shock[11]
        dist11 = people[:distance11]
      else
        val11 = -10^10
        dist11 = 999
      end
      choice = 0;
print("Entrant size:" , maximum(size(entrants)), "\n")
      if maximum(size(entrants)) <= 1 # there isn't a two-d array with size 0.
        print("In first loop part \n")
        chosen = indmax([val1.data val2.data val3.data val4.data val5.data val6.data val7.data val8.data val9.data val10.data val11.data]) # returns the *index* of the max element in the collection
      elseif maximum(size(entrants)) > 1 # at least one entrant
        cval = maximum([val1.data val2.data val3.data val4.data val5.data val6.data val7.data val8.data val9.data val10.data val11.data])
        mindist = minimum([dist1.data dist2.data dist3.data dist4.data dist5.data dist6.data dist7.data dist8.data dist9.data dist10.data dist11.data])
        chosen = indmax([val1.data val2.data val3.data val4.data val5.data val6.data val7.data val8.data val9.data val10.data val11.data])
        entfids = [entrants[x] for x in 1:ent_length:maximum(size(entrants))]
        #=
        Future Fix:
        Right now not correcting for the fact that the entrant might be the closest.
        Idea: indmin([distances]) - closest guy.  If [val1 val2 ... val11][indmin] == cval, then
        closest one is the maximizer.  Then cval = cval - closest_c.  Now do comparison.
        =#
        if entfids[1] != 0
          for k = 1:size(entfids)[1]
            efid = entfids[k]
            eind = findfirst(entrants, efid)
            # generate the distance x beds interaction
            dist = distance(people[:ZIP_LAT], people[:ZIP_LONG], entrants[eind+4], entrants[eind+5])
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
      else
          # undefined...
      end

      if chosen <= 11
        choicemade = convert(Int64, eval(parse(frname*"[:fid$chosen].data[1]")))
      else
        choicemade = convert(Int64, entfids[choice])
      end
    return choicemade
end



function ftest(people::DataFrame, demandmodelparameters::Array{Float64, 2})
  # Just returns the sum, compare to ftest2 which returns the distance
  shock = rand(d, 1)
  distance_c  = modelparameters[1]
  distsq_c  = modelparameters[2]
  neoint_c  = modelparameters[3]
  soloint_c  = modelparameters[4]
  closest_c  = modelparameters[5]
  distbed_c = modelparameters[6]

  return  people[:distance1]*distance_c + people[:distsq1]*distsq_c + people[:SoloIntermediate1]*soloint_c + people[:NeoIntensive1]*neoint_c + people[:closest1]*closest_c + people[:dist_bed1]*distbed_c + shock, people[:distance1]
end

function ftest2(people::DataFrame, demandmodelparameters::Array{Float64, 2})
  shock = rand(d, 1)
  distance_c  = modelparameters[1]
  distsq_c  = modelparameters[2]
  neoint_c  = modelparameters[3]
  soloint_c  = modelparameters[4]
  closest_c  = modelparameters[5]
  distbed_c = modelparameters[6]

  return  people[:distance1]*distance_c + people[:distsq1]*distsq_c + people[:SoloIntermediate1]*soloint_c + people[:NeoIntensive1]*neoint_c + people[:closest1]*closest_c + people[:dist_bed1]*distbed_c + shock
end

# x = [12, 17, 11, 5, 13, 16]
# convert(Vector, people[x])
mat = [ convert(Vector, people[12]) convert(Vector, people[17]) convert(Vector, people[11]) convert(Vector, people[5]) convert(Vector, people[13]) convert(Vector, people[16])]


function ftest3(people::Array{Float64,2}, demandmodelparameters::Array{Float64, 2})
  rows = size(people)[1]
  shock = rand(d, rows)
  return sum(repmat(demandmodelparameters, rows, 1).*people) +shock
end

function ftest4(people::Array{Float64,2}, demandmodelparameters::Array{Float64, 2})
  rows = size(people)[1]
  outp = zeros(rows)
  for i = 1:rows
    outp[i] = [demandmodelparameters*people[i,:]' + rand(d, 1)][1]
  end
  return outp
end

function ftest5(people::Array{Float64,2}, demandmodelparameters::Array{Float64, 2})
  rows = size(people)[1]
  outp = zeros(rows)
  for i = 1:rows
    outp[i] = [demandmodelparameters*people[i,:]'][1] + rand(d, 1)[1]
  end
  return outp
end

# Something like this should work....  That is, be faster...

function ftest6(inpmat::Array{Float64,2}, modelparameters::Array{Float64, 2})
  newmat = Array{Float64,1}(size(inpmat, 1))
  for i in 1:size(people,1)
      for j in 1:size(modelparameters,2)
          newmat[i] += inpmat[i,j]*modelparameters[j]
      end
      newmat[i]+= rand(d,1)[1]
  end
  return newmat
end

function ftest7(inpmat::Array{Float64,2}, modelparameters::Array{Float64, 2})
  newmat = Array{Float64,1}(size(inpmat, 1))
  rands = rand(d, size(inpmat,1))
  for i in 1:size(inpmat,1)
      for j in 1:size(modelparameters,2)
          newmat[i] += inpmat[i,j]*modelparameters[j]
      end
      newmat[i]+= rands[i]
  end
  return newmat
end

# Ok - this one is actually fast:
# 0.034627 seconds (15 allocations: 8.677 MB, 17.66% gc time)

function ftest8(inpmat::Array{Float64,2}, modelparameters::Array{Float64, 2})
  rands = rand(d, size(inpmat,1))
  return inpmat*modelparameters' + rands
end

# Now try to start with a dataframe and get the same.
# 0.041113 seconds (24 allocations: 34.707 MB, 17.63% gc time) - with ONE matrix, eg. mat1 only
# 0.634174 seconds (218 allocations: 381.779 MB, 50.13% gc time) - with 11 matrices
function imptest(peo::DataFrame, modelparameters::Array{Float64, 2}; siz = size(peo,1) , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176] )
  siz = size(people, 1)
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
#  rands = rand(d, size(people, 1), 11)
#  return mat1*modelparameters' + rands[:,1], mat2*modelparameters' + rands[:,2], mat3*modelparameters' + rands[:,3], mat4*modelparameters' + rands[:,4], mat5*modelparameters' + rands[:,5], mat6*modelparameters' + rands[:,6], mat7*modelparameters' + rands[:,7], mat8*modelparameters' + rands[:,8], mat9*modelparameters' + rands[:,9], mat10*modelparameters' + rands[:,10], mat11*modelparameters' + rands[:,11]
  return mat1 , mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9, mat10, mat11
end

# Little difference made by preallocating rands:
# 0.674994 seconds (246 allocations: 413.593 MB, 47.70% gc time)
# Effectively no difference made by returning [mat1 mat2 ... ] instead of "return mat1, mat2, ...":
# 0.684573 seconds (274 allocations: 445.408 MB, 48.64% gc time)
function imptest11(peo::DataFrame, modelparameters::Array{Float64, 2}; maxfid = 11, siz = size(peo,1) , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176] )
  siz = size(people, 1)
  rands = rand(d, size(people, 1), maxfid)
  mat1 = [ convert(Vector{Float64}, peo[ind[1]]) convert(Vector{Float64}, peo[ind[2]]) convert(Vector{Float64}, peo[ind[3]]) convert(Vector{Float64}, peo[ind[4]]) convert(Vector{Float64}, peo[ind[5]]) convert(Vector{Float64}, peo[ind[6]])]*modelparameters' + rands[:,1]
  mat2 = [ convert(Vector{Float64}, peo[iind[1]]) convert(Vector{Float64}, peo[iind[2]]) convert(Vector{Float64}, peo[iind[3]]) convert(Vector{Float64}, peo[iind[4]]) convert(Vector{Float64}, peo[iind[5]]) convert(Vector{Float64}, peo[iind[6]])]*modelparameters' + rands[:,2]
  mat3 = [ convert(Vector{Float64}, peo[iiind[1]]) convert(Vector{Float64}, peo[iiind[2]]) convert(Vector{Float64}, peo[iiind[3]]) convert(Vector{Float64}, peo[iiind[4]]) convert(Vector{Float64}, peo[iiind[5]]) convert(Vector{Float64}, peo[iiind[6]])]*modelparameters' + rands[:,3]
  mat4 = [ convert(Vector{Float64}, peo[ivnd[1]]) convert(Vector{Float64}, peo[ivnd[2]]) convert(Vector{Float64}, peo[ivnd[3]]) convert(Vector{Float64}, peo[ivnd[4]]) convert(Vector{Float64}, peo[ivnd[5]]) convert(Vector{Float64}, peo[ivnd[6]])]*modelparameters' + rands[:,4]
  mat5 = [ convert(Vector{Float64}, peo[vnd[1]]) convert(Vector{Float64}, peo[vnd[2]]) convert(Vector{Float64}, peo[vnd[3]]) convert(Vector{Float64}, peo[vnd[4]]) convert(Vector{Float64}, peo[vnd[5]]) convert(Vector{Float64}, peo[vnd[6]])]*modelparameters' + rands[:,5]
  mat6 = [ convert(Vector{Float64}, peo[vind[1]]) convert(Vector{Float64}, peo[vind[2]]) convert(Vector{Float64}, peo[vind[3]]) convert(Vector{Float64}, peo[vind[4]]) convert(Vector{Float64}, peo[vind[5]]) convert(Vector{Float64}, peo[vind[6]])]*modelparameters' + rands[:,6]
  mat7 = [ convert(Vector{Float64}, peo[viind[1]]) convert(Vector{Float64}, peo[viind[2]]) convert(Vector{Float64}, peo[viind[3]]) convert(Vector{Float64}, peo[viind[4]]) convert(Vector{Float64}, peo[viind[5]]) convert(Vector{Float64}, peo[viind[6]])]*modelparameters' + rands[:,7]
  mat8 = [ convert(Vector{Float64}, peo[viiind[1]]) convert(Vector{Float64}, peo[viiind[2]]) convert(Vector{Float64}, peo[viiind[3]]) convert(Vector{Float64}, peo[viiind[4]]) convert(Vector{Float64}, peo[viiind[5]]) convert(Vector{Float64}, peo[viiind[6]])]*modelparameters' + rands[:,8]
  mat9 = [ convert(Vector{Float64}, peo[ixnd[1]]) convert(Vector{Float64}, peo[ixnd[2]]) convert(Vector{Float64}, peo[ixnd[3]]) convert(Vector{Float64}, peo[ixnd[4]]) convert(Vector{Float64}, peo[ixnd[5]]) convert(Vector{Float64}, peo[ixnd[6]])]*modelparameters' + rands[:,9]
  mat10 = [ convert(Vector{Float64}, peo[xnd[1]]) convert(Vector{Float64}, peo[xnd[2]]) convert(Vector{Float64}, peo[xnd[3]]) convert(Vector{Float64}, peo[xnd[4]]) convert(Vector{Float64}, peo[xnd[5]]) convert(Vector{Float64}, peo[xnd[6]])]*modelparameters' + rands[:,10]
  mat11 = [ convert(Vector{Float64}, peo[xind[1]]) convert(Vector{Float64}, peo[xind[2]]) convert(Vector{Float64}, peo[xind[3]]) convert(Vector{Float64}, peo[xind[4]]) convert(Vector{Float64}, peo[xind[5]]) convert(Vector{Float64}, peo[xind[6]])]*modelparameters' + rands[:,11]
#  rands = rand(d, size(people, 1), 11)
#  return mat1*modelparameters' + rands[:,1], mat2*modelparameters' + rands[:,2], mat3*modelparameters' + rands[:,3], mat4*modelparameters' + rands[:,4], mat5*modelparameters' + rands[:,5], mat6*modelparameters' + rands[:,6], mat7*modelparameters' + rands[:,7], mat8*modelparameters' + rands[:,8], mat9*modelparameters' + rands[:,9], mat10*modelparameters' + rands[:,10], mat11*modelparameters' + rands[:,11]
  return [mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11]
end

# The result of:
# a1 = imptest11(people, modelparameters)
# @time mapslices(indmax, a1, 2)
# 0.448636 seconds (6.82 M allocations: 211.104 MB, 15.90% gc time)
# So this is the source of allocations later.

# Indices of the columns of interest for each fid:
#=
each in the order: [distance, distsq, SoloIntermediate, NeoIntensive, Closest, dist_bed]

for i in 1:11
    print("FID # ", i, "\n")
    print("[", people.colindex.lookup[eval(parse(":distance"*"$i"))], " ")
    print(people.colindex.lookup[eval(parse(":distsq"*"$i"))], " ")
    print(people.colindex.lookup[eval(parse(":SoloIntermediate"*"$i"))], " ")
    print(people.colindex.lookup[eval(parse(":NeoIntensive"*"$i"))], " ")
    print(people.colindex.lookup[eval(parse(":closest"*"$i"))], " ")
    print(people.colindex.lookup[eval(parse(":dist_bed"*"$i"))], "]", "\n ")
    print("********\n")
end
# List all fids
for i in 1:11
    print(people.colindex.lookup[eval(parse(":fid"*"$i"))], " ")
end
allfids = [2 18 34 50 66 82 98 114 130 146 162]

fid1: [12 17 11 5 13 16]
fid2: [28 33 27 21 29 32]
fid3: [44 49 43 37 45 48]
fid4: [60 65 59 53 61 64]
FID5: [76 81 75 69 77 80]
FID6: [92 97 91 85 93 96]
FID7: [108 113 107 101 109 112]
FID8: [124 129 123 117 125 128]
FID9: [140 145 139 133 141 144]
FID10: [156 161 155 149 157 160]
FID11: [172 177 171 165 173 176]

=#

# 0.602147 seconds (6.82 M allocations: 257.381 MB, 17.19% gc time)
function imptest2(peo::DataFrame, modelparameters::Array{Float64, 2}; indices = [12 17 11 5 13 16], iindices = [28 33 27 21 29 32])
  mat1 = [ convert(Vector{Float64}, peo[indices[1]]) convert(Vector{Float64}, peo[indices[2]]) convert(Vector{Float64}, peo[indices[3]]) convert(Vector{Float64}, peo[indices[4]]) convert(Vector{Float64}, peo[indices[5]]) convert(Vector{Float64}, peo[indices[6]])]*modelparameters' + rand(d, size(mat, 1))
  mat2 = [ convert(Vector{Float64}, peo[iindices[1]]) convert(Vector{Float64}, peo[iindices[2]]) convert(Vector{Float64}, peo[iindices[3]]) convert(Vector{Float64}, peo[iindices[4]]) convert(Vector{Float64}, peo[iindices[5]]) convert(Vector{Float64}, peo[iindices[6]])]*modelparameters' + rand(d, size(mat, 1))
  return mapslices(indmax, [mat1 mat2], 2)
end

# 0.469547 seconds (2.27 M allocations: 130.135 MB, 68.94% gc time)
function imptest3(peo::DataFrame, modelparameters::Array{Float64, 2}; indices = [12 17 11 5 13 16], iindices = [28 33 27 21 29 32])
  outp = zeros(size(peo,1))
  mat1 = [ convert(Vector{Float64}, peo[indices[1]]) convert(Vector{Float64}, peo[indices[2]]) convert(Vector{Float64}, peo[indices[3]]) convert(Vector{Float64}, peo[indices[4]]) convert(Vector{Float64}, peo[indices[5]]) convert(Vector{Float64}, peo[indices[6]])]*modelparameters' + rand(d, size(mat, 1))
  mat2 = [ convert(Vector{Float64}, peo[iindices[1]]) convert(Vector{Float64}, peo[iindices[2]]) convert(Vector{Float64}, peo[iindices[3]]) convert(Vector{Float64}, peo[iindices[4]]) convert(Vector{Float64}, peo[iindices[5]]) convert(Vector{Float64}, peo[iindices[6]])]*modelparameters' + rand(d, size(mat, 1))
  for i = 1:size(peo, 1)
    outp[i] = indmax([mat1[i], mat2[i]])
  end
  return outp
end

# 0.269869 seconds (3.03 M allocations: 182.188 MB, 28.28% gc time)
function imptest31(peo::DataFrame, modelparameters::Array{Float64, 2}; indices = [12 17 11 5 13 16], iindices = [28 33 27 21 29 32], iiindices = [44 49 43 37 45 48])
  outp = zeros(size(peo,1))
  mat1 = [ convert(Vector{Float64}, peo[indices[1]]) convert(Vector{Float64}, peo[indices[2]]) convert(Vector{Float64}, peo[indices[3]]) convert(Vector{Float64}, peo[indices[4]]) convert(Vector{Float64}, peo[indices[5]]) convert(Vector{Float64}, peo[indices[6]])]*modelparameters' + rand(d, size(mat, 1))
  mat2 = [ convert(Vector{Float64}, peo[iindices[1]]) convert(Vector{Float64}, peo[iindices[2]]) convert(Vector{Float64}, peo[iindices[3]]) convert(Vector{Float64}, peo[iindices[4]]) convert(Vector{Float64}, peo[iindices[5]]) convert(Vector{Float64}, peo[iindices[6]])]*modelparameters' + rand(d, size(mat, 1))
  mat3 = [ convert(Vector{Float64}, peo[iiindices[1]]) convert(Vector{Float64}, peo[iiindices[2]]) convert(Vector{Float64}, peo[iiindices[3]]) convert(Vector{Float64}, peo[iiindices[4]]) convert(Vector{Float64}, peo[iiindices[5]]) convert(Vector{Float64}, peo[iiindices[6]])]*modelparameters' + rand(d, size(mat, 1))
  for i = 1:size(peo, 1)
    outp[i] = indmax([mat1[i], mat2[i], mat3[i]])
  end
  return countmap(outp)
end

# test results: 0.942355 seconds (9.09 M allocations: 575.473 MB, 40.63% gc time)
# sometimes it takes slightly less time but average over 10 runs probably about 0.81 seconds
function imptest32(peo::DataFrame, modelparameters::Array{Float64, 2}; siz = size(peo,1) , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176] )
  outp = zeros(siz)
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
  for i = 1:size(peo, 1)
    outp[i] = indmax([mat1[i], mat2[i], mat3[i], mat4[i], mat5[i], mat6[i], mat7[i], mat8[i], mat9[i], mat10[i], mat11[i]])
  end
  return countmap(outp)
end


# Now the challenge is to add distance and permit entry

###

function imptest33(peo::DataFrame, modelparameters::Array{Float64, 2}; siz = size(peo,1) , ind = [12 17 11 5 13 16], iind = [28 33 27 21 29 32], iiind = [44 49 43 37 45 48], ivnd = [60 65 59 53 61 64], vnd = [76 81 75 69 77 80], vind = [92 97 91 85 93 96], viind = [108 113 107 101 109 112], viiind = [124 129 123 117 125 128], ixnd = [140 145 139 133 141 144], xnd = [156 161 155 149 157 160], xind = [172 177 171 165 173 176], fidnd = [2 18 34 50 66 82 98 114 130 146 162] )
  outp = zeros(siz)
  allfids = [convert(Vector{Int64}, peo[fidnd[1]]) convert(Vector{Int64}, peo[fidnd[2]]) convert(Vector{Int64}, peo[fidnd[3]]) convert(Vector{Int64}, peo[fidnd[4]]) convert(Vector{Int64}, peo[fidnd[5]]) convert(Vector{Int64}, peo[fidnd[6]]) convert(Vector{Int64}, peo[fidnd[7]]) convert(Vector{Int64}, peo[fidnd[8]]) convert(Vector{Int64}, peo[fidnd[9]]) convert(Vector{Int64}, peo[fidnd[10]]) convert(Vector{Int64}, peo[fidnd[11]])]
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
  for i = 1:size(peo, 1) #change the version above: list all fids for each person, then take i-th (person) row, indmax columns
    outp[i] = allfids[i, indmax([mat1[i], mat2[i], mat3[i], mat4[i], mat5[i], mat6[i], mat7[i], mat8[i], mat9[i], mat10[i], mat11[i]])]
  end
  return countmap(outp)
end
