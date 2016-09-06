"""
correct approach - these are going to be probabilities, not realized choices.
Change the demand model thing - make two functions: compute deterministic component
of utility, then choice probs as a function of that, but use same output as input
for adding random shocks and perturbing to figure out actual choice - will cut allocations.
This will be easier than expected.

"""

# choicerow/people  has the form: [identity] ∪ [fid, NeoIntensive, Solo Intermediate, distance, Is Closest?, Selected?, distance × bed, distance², amount charged] × (# facilities) ⋃ [Patient Zip, DRG, medicaid, Private insurance, Zip Lat, Zip Long]

# Takes about 0.09 - 0.12 seconds per call: 0.124532 seconds (334 allocations: 82.751 MB, 37.50% gc time)

function DetUtil(peo::Matrix, modelparameters::Array{Float64, 2}; ziploc = 101, drgloc = 102,  fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92] , fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24], ind = [5 9 3 4 6 8 ], iind = [14 18 12 13 15 17 ], iiind = [23 27 21 22 24 26 ], ivnd = [32 36 30 31 33 35 ], vnd = [41 45 39 40 42 44 ], vind = [50 54 48 49 51 53 ], viind = [59 63 57 58 60 62 ], viiind = [68 72 66 67 69 71 ], ixnd = [77 81 75 76 78 80 ], xnd = [86 90 84 85 87 89 ], xind = [95 99 93 94 96 98 ] )
# Computed utilities
# Use as input to WTP as market shares and to Demand Model by adding random shocks
# This is inefficient because it does it for *everyone*, which we really don't need.
  mat1 = peo[:,ind[1:6]]*modelparameters'
  mat2 = peo[:,iind[1:6]]*modelparameters'
  mat3 = peo[:,iiind[1:6]]*modelparameters'
  mat4 = peo[:,ivnd[1:6]]*modelparameters'
  mat5 = peo[:,vnd[1:6]]*modelparameters'
  mat6 = peo[:,vind[1:6]]*modelparameters'
  mat7 = peo[:,viind[1:6]]*modelparameters'
  mat8 = peo[:,viiind[1:6]]*modelparameters'
  mat9 = peo[:,ixnd[1:6]]*modelparameters'
  mat10 = peo[:,xnd[1:6]]*modelparameters'
  mat11 = peo[:,xind[1:6]]*modelparameters'
  return hcat( peo[:, ziploc], peo[:, drgloc], mat1, peo[:, fidnd[1]], mat2, peo[:, fidnd[2]], mat3, peo[:, fidnd[3]], mat4, peo[:, fidnd[4]], mat5, peo[:, fidnd[5]], mat6, peo[:, fidnd[6]], mat7, peo[:, fidnd[7]],  mat8, peo[:, fidnd[8]],  mat9, peo[:, fidnd[9]], mat10, peo[:, fidnd[10]], mat11, peo[:, fidnd[11]])
end

# I DON'T want to call DetUtil inside this function - the whole point of DetUtil is to avoid Multiple calls to that.
# test = DetUtil(pinsured, privatedemandmodelparameters)

# 0.221892 seconds (5.12 M allocations: 263.652 MB, 20.08% gc time)
function ComputeWTP(utils::Matrix ; ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
  output = zeros(size(utils, 1), size(fidlocs,1) )
  for i = 1:size(utils,1) # over rows of return type
    utils[i, ulocs[:]] = exp(utils[i, ulocs[:]])
    interim = 0
    for j in fidlocs
      if utils[i,j] != 0 # the facility is not absent
#        utils[i, j-1] = exp(utils[i,j-1])
        interim += utils[i,j-1] # add the previous value
      else
        utils[i, j-1] = 0
      end
    end
    utils[i,ulocs[:]] = utils[i, ulocs[:]]./interim
  end
  return utils
end # of ComputeWTP

# This is a tester - does the function ComputeWTP return a set of probabilities or not?
# It does to that, to within floating point errors, basically.
function CheckWTP(peo::Matrix; params = privatedemandmodelparameters, ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23])
  util = DetUtil(peo, params)
  probs = ComputeWTP(util)
  checkvec = sum( probs[:, ulocs[:]], 2)
  return checkvec[(checkvec.>1.0 + eps())|(checkvec.<1.0-eps())]
end

# Call this after DetUtil and ComputeWTP
# test = DetUtil(pinsured, privatedemandmodelparameters)
# tst_probs = ComputeWTP(test)
# Output of ComputeWTP is (Zip, DRG, ) ∪ (Utility, Hospital) × 12
# DRGs = 385 386 387 388 389 390 391
# IF instead I take JUST the unique zip code DRG sets, then I can get the time down from 7 seconds to 0.3-0.4 seconds.
function MapWTP(comp_wtp::Matrix ; pziploc = 1, pdrgloc = 2, zipcodes = TXzips, fids = allfids, drg = DRGs, modparams = privatedemandmodelparameters, ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
  # Creates the Output:
  # This block takes kind of a long time: 0.028215 seconds (132.51 k allocations: 39.350 MB, 18.91% gc time)
  output = zeros( size(zipcodes,1)*size(drg,2)+1, size(fids, 1)+3)
  for i = 1:size(fids, 1)
    output[1, i+2] = fids[i]
  end
  for i = 1:size(zipcodes, 1)
    for j =1:size(drg, 2)
      output[7*(i-1)+j+1, 1] = zipcodes[i] # first column is zip
      output[7*(i-1)+j+1, 2] = drg[j] # second column is drg.
      # third column reports whether found or not.
    end
  end
  # Map the comp_wtp to the matrix of values
  for i = 1:size(comp_wtp,1)
    st = findfirst(output[:,1], comp_wtp[i, pziploc]) # find first occurrence of zipcode
    for j = 0:size(drg, 2)-1 # this is inefficient - I know the order and the number of DRGs every time.  Maybe hardcode that?
      if output[st+j, 2] == comp_wtp[i, pdrgloc]
        if output[st+j, 3] == 0 # not found before
          for k in fidlocs # indexes columns containing fids
            hos = findfirst(output[1,3:end], comp_wtp[i,k]) #find the fid in the first row of output
            if hos > 0
              if output[st+j,k] == 0
                output[st+j, k] += comp_wtp[i, k-1]
              else
                println("Row ", i, " column ", k, " problem ")
              end
            end
          end
        end
        output[st+j, 3] += 1
      end
    end
  end
  # This is sure to be horribly slow
  return output
end # of WTP



function MapWTP2(comp_wtp::Matrix ; pziploc = 1, pdrgloc = 2, zipcodes = TXzips, fids = allfids, drg = DRGs, ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
  # Creates the Output:
  # This block takes kind of a long time: 0.028215 seconds (132.51 k allocations: 39.350 MB, 18.91% gc time)
  output = zeros( size(zipcodes,1)*size(drg,2)+1, size(fids, 1)+3)
  for i = 1:size(fids, 1)
    output[1, i+3] = fids[i]
  end
  for i = 1:size(zipcodes, 1)
    for j =1:size(drg, 2)
      output[7*(i-1)+j+1, 1] = zipcodes[i] # first column is zip
      output[7*(i-1)+j+1, 2] = drg[j] # second column is drg.
      # third column reports whether found or not.
    end
  end
  # Map the comp_wtp to the matrix of values
  curr = 1
  for i = 1:size(comp_wtp,1)
    st = findfirst(output[curr:end,1], comp_wtp[i, pziploc]) # find first occurrence of zipcode
    j = convert(Int, comp_wtp[i, pdrgloc]-384)
    # if j > 7 # can't actually be greater than 7 by construction in stata
    #   println("Messed up DRG code row ", i, " ", comp_wtp[i,j])
    #   j = 0
    # end
    if output[st+j, 3] == 0 # not found before
      for k in fidlocs # indexes columns containing fids
        hos = findfirst(output[1,3:end], comp_wtp[i,k]) #find the fid in the first row of output
        if hos > 0
          output[st+j, k] += comp_wtp[i, k-1]
        end
      end
    end
    output[st+j, 3] += 1
    curr = st
  end
  return output
end # of MapWTP2



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
        for k = 3:size(output, 2) # indexes columns / hospitals
          if output[1, k] == a1[i, pfidloc] # this is wrong.
            output[st+j, k] += 1
          end
        end
      end
    end
  end
end




output = zeros( size(TXzips,1)*size(DRGs,2)+1, size(allfids, 1)+2)
for i = 1:size(allfids, 1)
  output[1, i+2] = allfids[i]
end
for i = 1:size(TXzips, 1)
  for j =1:size(DRGs, 2)
    output[7*(i-1)+j+1, 1] = TXzips[i]
    output[7*(i-1)+j+1, 2] = DRGs[j]
  end
end
