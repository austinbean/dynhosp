# Use the parameters of the demand model to estimate choices:

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
# sample entrants3 = [99999 1 0 120 32.96  -96.8385 888888 0 1 120 32.96  -96.8385 77777 0 0 120 31.96  -97.8385]
# with 3 entrants this takes 0.05 seconds.



function EntrantsU(peo::Matrix, entrants::Array{Float64, 2}, modelparameters::Array{Float64, 2};
                   dist_μ = 0,
                   dist_σ = 1,
                   dist_ξ = 0,
                   d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ),
                   persloc=[25,26],
                   entsize = 6,
                   entnum = convert(Int, size(entrants, 2)/entsize))
  siz = size(peo,1)
  rands = rand(d, siz, entnum)
  entvals = zeros(siz, entnum)
  plocs = peo[:,persloc[:]] # this format is: LATITUDE, LONGITUDE
  entfids = entrants[:, 1:entsize:end] # selects all fids.
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
  outp = hcat(utils, map(x->entfids[x], ind2sub((siz, size(entfids,1)), vec(inds))[2]))
  return outp  # note that due to the randomization, this will generally not return -999, but -999 + rand
end


# Call DetUtil first, then this.
function DemandModel(detutil::Matrix, modelparameters::Array{Float64, 2}, entrants::Array{Float64, 2};
                      dist_μ = 0,
                      dist_σ = 1,
                      dist_ξ = 0,
                      d = Distributions.GeneralizedExtremeValue(dist_μ, dist_σ, dist_ξ),
                      ziploc = 1,
                      drgloc = 2,
                      entsize = 6,
                      entnum = convert(Int, size(entrants, 2)/entsize),
                      siz = size(detutil,1),
                      fidnd = [2; 11; 20; 29; 38; 47; 56; 65; 74; 83; 92],
                      ulocs = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23],
                      fidlocs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
# Computed utilities + error
  rand_el = Array{Float64}(siz, 11)
  if size(entrants, 2) > 1
    entutil = EntrantsU(detutil, entrants, modelparameters)
    both = hcat(detutil[:,ulocs[:]] + rand!(d, rand_el), entutil[:,1])
    fids = hcat(detutil[:, fidlocs[:]], entutil[:,2])
    vals, inds = findmax( both , 2)
    # I need to recover the fid of the max of [detutil entutil by row ]
    outp = map((i,x)->fids[i,x], 1:siz, ind2sub((size(detutil,1),11 + 1), vec(inds) )[2] )
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
















#end
