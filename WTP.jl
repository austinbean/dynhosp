
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

function WTP(peo::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2}; entsize = 6, entnum = convert(Int, size(entrants, 2)/entsize), siz = size(peo,1), pziploc = 1, pdrgloc = 2, pfidloc = 3)
  for i = 1:size(peo,1)
    st = findfirst(output[:,1], peo[i, pziploc]) # find first occurrence of zipcode
    matched = false
    for j = 0:size(DRGs, 2)-1
      if output[st+j, k] == peo[i, pdrgloc]
        for k = 1:size(output, 2) # indexes columns
          if output[st+j, k] == peo[i, pfidloc]
            output[st+j, k] += 1
            matched = true
            # I want to break out of this loop after matching
          end
        end
      end
    end
  end
  return output
end # of WTP




pziploc = 1
pdrgloc = 2
pfidloc = 3
for i = 1:size(a1,1)
  st = findfirst(output[:,1], a1[i, pziploc]) # find first occurrence of zipcode
  if st == 0
    println(a1[i,pziploc])
  else
    for j = 0:size(DRGs, 2)-1
      if output[st+j, 2] == a1[i, pdrgloc]
        for k = 1:size(output, 2) # indexes columns
          if output[st+j, k] == a1[i, pfidloc]
            output[st+j, k] += 1
          end
        end
      end
    end
  end
end



# THIS IS A PROBLEM - which hospital fids?  How will that work? It has to be all of them

output = zeros( size(TXzips,1)*size(DRGs,2)+1, size(allfids, 1)+2)
for i = 1:size(allfids, 1)
  output[1, i+2] = allfids[i]
end
for i = 1:size(TXzips, 1)
  for j =1:size(DRGs, 2)
    println(7*(i-1)+j+1 )
    output[7*(i-1)+j+1, 1] = TXzips[i]
    output[7*(i-1)+j+1, 2] = DRGs[j]
  end
end
