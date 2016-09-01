
function EntrantsWTP(peo::Matrix, entrants::Array{Float64, 2}, modelparameters::Array{Float64, 2};  persloc = [ 105 106], entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize))
  siz = size(peo,1)
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

function WTP(peo::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92] , ind = [5 9 3 4 6 8 ], iind = [14 18 12 13 15 17 ], iiind = [23 27 21 22 24 26 ], ivnd = [32 36 30 31 33 35 ], vnd = [41 45 39 40 42 44 ], vind = [50 54 48 49 51 53 ], viind = [59 63 57 58 60 62 ], viiind = [68 72 66 67 69 71 ], ixnd = [77 81 75 76 78 80 ], xnd = [86 90 84 85 87 89 ], xind = [95 99 93 94 96 98 ] )


  if size(entrants, 2) > 1
    entfids = convert(Vector{Int64}, [entrants[x] for x in 1:entsize:size(entrants,2)])'
    allfids = [peo[:,fidnd[1:11]] repmat(entfids, siz, 1)] #maybe fill(entfids, siz) would work faster?  Maybe speed is the same but allocations lower.
    entutil = EntrantsU(peo, entrants, modelparameters)
    vals, inds = findmax([mat1 mat2 mat3 mat4 mat5 mat6 mat7 mat8 mat9 mat10 mat11 entutil[:,1]] , 2)
    # Fuck-up is right here.
    outp = map( (i,x)->allfids[i,x], collect(1:size(mat1,1)), ind2sub((size(mat1,1),12), vec(inds) )[2] )
  else #  no entrants


  end
return outp
end
