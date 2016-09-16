"""
correct approach - these are going to be probabilities, not realized choices.
Change the demand model thing - make two functions: compute deterministic component
of utility, then choice probs as a function of that, but use same output as input
for adding random shocks and perturbing to figure out actual choice - will cut allocations.
This will be easier than expected.

"""


# choicerow/people  has the form: [identity] ∪ [fid, NeoIntensive, Solo Intermediate, distance, Is Closest?, Selected?, distance × bed, distance², amount charged] × (# facilities) ⋃ [Patient Zip, DRG, medicaid, Private insurance, Zip Lat, Zip Long]

# Takes about 0.09 - 0.12 seconds per call: 0.124532 seconds (334 allocations: 82.751 MB, 37.50% gc time)

function FindCorrect(peo::Matrix; ziploc = 101, drgloc = 104, persloc = [ 107 108] )
  checklocs = true
  try
     peo[1,ziploc]
   catch errz
     if isa(errz, BoundsError)
       checklocs = false
       println("Size Changed")
       println("ZIPLOC PROBABLY HERE: ")
       println(findfirst(peo[1,:], 75001))
       println("*******************")
     elseif peo[1,ziploc] != 75001
      println("FIX THE ZIP CODE LOCATION IN DETUTIL")
      println("IT'S PROBABLY HERE: ")
      println(findfirst(peo[1,:], 75001))
      println("*******************")
      checklocs = false
    end
  end
    try
      peo[1,drgloc]
    catch errz
      if isa(errz, BoundsError)
        println("Size Changed")
        println("DRGLOC PROBABLY HERE: ")
        println(findfirst(peo[1,:], 391))
        println("*******************")
        checklocs = false
      elseif peo[1, drgloc ]!= 391
        println("FIX THE DRG CODE LOCATION IN DETUTIL")
        println("DRGLOC PROBABLY HERE: ")
        println(findfirst(peo[1,:], 391))
        println("*******************")
        checklocs = false
      end
    end
    try
      peo[1,persloc[:]]
    catch errz
      if isa(errz, BoundsError)
        println("size Changed")
        println("LOCATION PROBABLY HERE: ")
        println(findfirst(peo[1,:], convert(Float32, 32.960049)))
        println(findfirst(peo[1,:], convert(Float32, -96.838524)))
        println("*******************")
        checklocs = false
      elseif peo[1,persloc[:]] != [convert(Float32,32.960049)	convert(Float32, -96.838524)]
        println("FIX THE LAT/LONG CODE LOCATION IN DETUTIL")
        println("IT'S PROBABLY HERE: ")
        println(findfirst(peo[1,:], convert(Float32, 32.960049)))
        println(findfirst(peo[1,:], convert(Float32, -96.838524)))
        println("*******************")
        checklocs = false
      end
    end
    if checklocs
      println("Locations Correct")
    end
end




function DetUtil(peo::Matrix, modelparameters::Array{Float64, 2}; ziploc = 101, drgloc = 104, persloc = [107 108],  fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92] , fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24], ind = [5 9 3 4 6 8 ], iind = [14 18 12 13 15 17 ], iiind = [23 27 21 22 24 26 ], ivnd = [32 36 30 31 33 35 ], vnd = [41 45 39 40 42 44 ], vind = [50 54 48 49 51 53 ], viind = [59 63 57 58 60 62 ], viiind = [68 72 66 67 69 71 ], ixnd = [77 81 75 76 78 80 ], xnd = [86 90 84 85 87 89 ], xind = [95 99 93 94 96 98 ] )
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
  vals = hcat( peo[:, ziploc], peo[:, drgloc], mat1, peo[:, fidnd[1]], mat2, peo[:, fidnd[2]], mat3, peo[:, fidnd[3]], mat4, peo[:, fidnd[4]], mat5, peo[:, fidnd[5]], mat6, peo[:, fidnd[6]], mat7, peo[:, fidnd[7]],  mat8, peo[:, fidnd[8]],  mat9, peo[:, fidnd[9]], mat10, peo[:, fidnd[10]], mat11, peo[:, fidnd[11]], peo[:, persloc[:]])
  # Adding the next loop doubles - triples the time per function evaluation: 0.08 seconds to 0.24 seconds for ≈ 150,000 rows.
  for i = 1:size(vals,1)
    for j in fidlocs
      if vals[i,j] == 0
        vals[i, j-1] = -99 # prevent missing facilities from being chosen by assigning them large negative utility.
      end
    end
  end
  return vals
end




# 0.221892 seconds (5.12 M allocations: 263.652 MB, 20.08% gc time)
function ComputeWTP(utils::Matrix ; ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
  output = zeros(size(utils))
  for i = 1:size(utils,1) # over rows of return type
    interim = 0
    for j in 1:size(utils,2)
      if (j == 1)|(j == 2)
        output[i,j] = utils[i,j] # be careful - I'll be this *doesn't* allocate a new array, just a set of pointers or something.
      elseif j in fidlocs
        output[i,j] = utils[i, j]
        if utils[i,j] != 0 # the facility is not absent
          output[i,j-1] = exp(utils[i,j-1])
          interim += exp(utils[i,j-1]) # add the previous value
        else
          output[i, j-1] = 0 # be careful - I'll be this *doesn't* allocate a new array
        end
      end
    end
    output[i,ulocs[:]] = output[i, ulocs[:]]./interim
  end
  return output
end # of ComputeWTP


# This is a tester - does the function ComputeWTP return a set of probabilities or not?
# It does to that, to within floating point errors, basically.
function CheckWTP(peo::Matrix; params = privatedemandmodelparameters, ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23])
  util = DetUtil(peo, params)
  probs = ComputeWTP(util)
  checkvec = sum( probs[:, ulocs[:]], 2)
  return checkvec[(checkvec.>1.0 + eps())|(checkvec.<1.0-eps())]
end

#testprobs = CheckWTP(pinsured)

# Output of ComputeWTP is (Zip, DRG, ) ∪ (Utility, Hospital) × 12
# DRGs = 385 386 387 388 389 390 391
# IF instead I take JUST the unique zip code DRG sets, then I can get the time down from 7 seconds to 0.3-0.4 seconds.


# This function tests whether some of the computed probs in ComputeWTP are 1 - this will create a problem in MapWTP as the valuation will be ∞

function ChoiceCount(comp_wtp::Matrix; fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
  ones_f = Array{Float64,1}()
  onechoice = Array{Float64,1}()
  nochoice = Array{Float64,1}()
  for i = 1:size(comp_wtp,1)
      counting = 0
      for j in fidlocs
        if comp_wtp[i,j] != 0
          counting += 1
        end
        if (comp_wtp[i,j-1] >= 1 - eps())&&(comp_wtp[i,j-1] <= 1+eps())
          println("Found a 1 ", i, "  ", j-1)
          push!(ones_f, i)
        end
      end
      if counting == 1
        println("One Choice at ", i)
        push!(onechoice, i)
      elseif counting == 0
        println("Zero Choices at ", i)
        push!(nochoice, i)
      end
  end
  return onechoice, nochoice, ones_f
end

#onecn, noch, ones_f = ChoiceCount(testWTP)

TXzps = DataFrames.readtable(pathprograms*"TXzipsonly.csv", header = false)
TXzips = convert(Vector, TXzps[:x1])
txfd = DataFrames.readtable(pathprograms*"TXfidsonly.csv", header = true)
allfids = convert(Vector, txfd[:fid])


# Record format:
# :fid1
 :NeoIntensive1
 :SoloIntermediate1
 :distance1
 :hosplat1
 :hosplong1
 :distsq1
 :closest1
 :distbed1
 privatedemandmodelparameters = [privatedistance_c privatedistsq_c privateneoint_c privatesoloint_c privateclosest_c privatedistbed_c]
ch1 =  [6 9 4 5 10 11]
ch2 = [15 18 13 14 19 20]
ch3 = [24 27 22 23 28 29]
ch4 = [33 36 31 32 37 38]
ch5 = [42 45 40 41 46 47]
ch6 = [51 54 49 50 55 56]
ch7 = [60 63 58 59 64 65]
ch8 = [69 72 67 68 73 74]
ch9 = [78 81 76 77 82 83]
ch10 = [87 90 85 86 91 92]

function EasyWTP(zipdrg::Matrix, modelparameters::Array{Float64,2};
                 fids = allfids,
                 drg = [385 386 387 388 389 390 391],
                 ch1 =  [6 9 4 5 10 11],
                 ch2 = [15 18 13 14 19 20],
                 ch3 = [24 27 22 23 28 29],
                 ch4 = [33 36 31 32 37 38],
                 ch5 = [42 45 40 41 46 47],
                 ch6 = [51 54 49 50 55 56],
                 ch7 = [60 63 58 59 64 65],
                 ch8 = [69 72 67 68 73 74],
                 ch9 = [78 81 76 77 82 83],
                 ch10 = [87 90 85 86 91 92])
  wtp1 = zipdrg[:, ch1[:]]*modelparameters'
  wtp2 = zipdrg[:, ch2[:]]*modelparameters'
  wtp3 = zipdrg[:, ch3[:]]*modelparameters'
  wtp4 = zipdrg[:, ch4[:]]*modelparameters'
  wtp5 = zipdrg[:, ch5[:]]*modelparameters'
  wtp6 = zipdrg[:, ch6[:]]*modelparameters'
  wtp7 = zipdrg[:, ch7[:]]*modelparameters'
  wtp8 = zipdrg[:, ch8[:]]*modelparameters'
  wtp9 = zipdrg[:, ch9[:]]*modelparameters'
  wtp10 = zipdrg[:, ch10[:]]*modelparameters'

end




# TODO: this is not doing exactly what I want it to, I think.  Check TX WTP Creator.do for discussion.  But there is reason to be concerned.
# The issue is that this sums over all of the WTP's, which is not what I want - or maybe it is but I want to be able to treat them separately at first.
# What do I want to do?  Do this by DRG.  That should remain constant over time.  But I only need one per type



function MapWTP(comp_wtp::Matrix ; pziploc = 1, pdrgloc = 2, zipcodes = TXzips, fids = allfids, drg = [385 386 387 388 389 390 391], ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
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
  for i = 1:size(comp_wtp,1)
    st = findfirst(output[:,1], comp_wtp[i, pziploc]) # find first occurrence of zipcode
    j = convert(Int, comp_wtp[i, pdrgloc]-385) # this is always the same
    if output[st+j, 3] == 0 # not found before
      for k in fidlocs # indexes columns containing fids
        hos = findfirst(output[1,3:end], comp_wtp[i,k]) #find the fid in the first row of output
        if hos > 0
          output[st+j, k] += log(1/(1-comp_wtp[i, k-1]))
        end
      end
    end
    output[st+j, 3] += 1
  end
  return output
end # of MapWTP2



function ReturnWTP(mapped_wtp::Matrix)
  return vcat( mapped_wtp[1,4:end], sum( mapped_wtp[2:end, 4:end].*mapped_wtp[2:end,3] ,1))
end


function CheckMaxWTP(wtp::Matrix)
  for i = 1:size(wtp,2)
    if wtp[2,i] == Inf
      println("Inf Found ", wtp[1,i])
      println(wtp[2,i])
    end
  end
  return maximum(wtp[2,:]), minimum(wtp[2,:])
end

#extremawtp = CheckMaxWTP(sumWTP)


WTPvals = DataFrames.readtable(pathprograms*"WTPTemplate.csv", header = true)
WTPvals = convert(Matrix, WTPvals)



#=
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
=#


# I need a WTP test to see which one of the two is working::
# choicerow/people  has the form: [identity] ∪ [fid, NeoIntensive, Solo Intermediate, distance, Is Closest?, Selected?, distance × bed, distance², amount charged] × (# facilities) ⋃ [Patient Zip, DRG, medicaid, Private insurance, Zip Lat, Zip Long]
# model parameters:      [distance distsq      neoint   soloint     closest distbed]


# Test the Deterministic Components:
#=
testmodelparameters = [-0.172602  0.00456541  0.962496  0.674374  0.53948  0.0171834]
#            [ID   FID   NEO SOLO DIS CLO SEL D×B D^2 $$$]
testpeople1 = [1.0 12719.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               2.0 12719.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               3.0 12719.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               4.0 12719.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               5.0 12719.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               6.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               7.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               8.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               9.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
              10.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524]

a1 = DetUtil(testpeople1, privatedemandmodelparameters)


# Test the Computation of WTP

testpeople2 = [1.0 12719.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               2.0 12719.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               3.0 12719.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               4.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               5.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               6.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               7.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               8.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               9.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
              10.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524]

testparameters2 = [1.0 2 1 1 1 1]
a2 = DetUtil(testpeople2, testparameters2)
b2 = ComputeWTP(a2)

# First row should be 1/2 1/2
# Second should be 1/3 1/3 1/3
# Third should be 1/12
# fourth should be NaN.


# Test the mapping of WTP
# This tested with WTP mapped just as prob directly, NOT as ln(1/1-prob)

testpeople3 = [1.0 12719.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 00000.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 000000.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 000000.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 00000.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               2.0 12719.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 000000.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.00000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 000000.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0000000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 00000.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               3.0 12719.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 391.0 0.0 1.0 32.96005 -96.838524;
               4.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 390.0 0.0 1.0 32.96005 -96.838524;
               5.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 390.0 0.0 1.0 32.96005 -96.838524;
               6.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 386.0 0.0 1.0 32.96005 -96.838524;
               7.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 386.0 0.0 1.0 32.96005 -96.838524;
               8.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 386.0 0.0 1.0 32.96005 -96.838524;
               9.0 12719.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 386.0 0.0 1.0 32.96005 -96.838524;
              10.0 12719.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 16122.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 30105.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.216088e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 856364.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136268e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.136007e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13105e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.135009e6 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 855094.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.13106e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1415013 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 75001.0 386.0 0.0 1.0 32.96005 -96.838524]

testparameters2 = [1.0 2 1 1 1 1]
a3 = DetUtil(testpeople3, testparameters2)
b3 = ComputeWTP(a3)
c3 = MapWTP(b3)
# sum(c4[:,3])==size(testpeople3, 1) # -> must return true.

=#



#=
# Slightly slower version of MapWTP

function MapWTP2(comp_wtp::Matrix ; pziploc = 1, pdrgloc = 2, zipcodes = TXzips, fids = allfids, drg = DRGs, modparams = privatedemandmodelparameters, ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23], fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
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


=#
