# reboot:
# include("/Users/austinbean/Desktop/dynhosp/Reboot.jl")
# include("/home/ubuntu/dynhosp/Reboot.jl")
# include("/home1/04179/abean/dynhosp/Reboot.jl")


#lis = addprocs(2)
#lis = addprocs() # for the 32 core Amazon machine.



push!(LOAD_PATH, "/Users/austinbean/Desktop/dynhosp")
push!(LOAD_PATH, "/dynhosp/dynhosp")
push!(LOAD_PATH, "/home/ubuntu/dynhosp/")
push!(LOAD_PATH, "/home1/04179/abean/dynhosp")


using ProjectModule
using DataFrames
using Distributions


#@everywhere # it seems that the data objects do *not* need to be defined everywhere.  That's weird.
begin
  # Figure out which machine I'm on
  dir = pwd()
  global pathdata = "";global pathpeople = "";global  pathprograms = "";
  if dir == "/Users/austinbean/Desktop/dynhosp"
    global pathdata = "/Users/austinbean/Google Drive/Annual Surveys of Hospitals/"
    global pathpeople = "/Users/austinbean/Google Drive/Texas Inpatient Discharge/"
    global pathprograms = "/Users/austinbean/Desktop/dynhosp/"
  elseif dir == "/dynhosp/dynhosp"
    global pathdata = dir*"/"
    global pathpeople = dir*"/"
    global pathprograms = dir*"/"
  elseif (dir == "/home/ubuntu/Notebooks") | (dir == "/home/ubuntu/dynhosp")
    global pathdata = "/home/ubuntu/dynhosp/"
    global pathpeople = "/home/ubuntu/dynhosp/"
    global pathprograms = "/home/ubuntu/dynhosp/"
  elseif (dir == "/home1/04179/abean/dynhosp")
    global pathdata = "/home1/04179/abean/dynhosp"
    global pathpeople = "/home1/04179/abean/dynhosp"
    global pathprograms = "/home1/04179/abean/dynhosp"
  end
  # Import the hospital data and convert to a matrix -
  println("Importing Hosp Data")
    dataf = DataFrames.readtable(pathdata*"TX Transition Probabilities.csv", header = true);
    for i in names(dataf)
        if ( typeof(dataf[i]) == DataArrays.DataArray{Float64,1} )
            dataf[DataFrames.isna(dataf[i]), i] = 0
        elseif (typeof(dataf[i]) == DataArrays.DataArray{Int64,1})
            dataf[DataFrames.isna(dataf[i]), i] = 0
        elseif typeof(dataf[i]) == DataArrays.DataArray{ByteString,1}
            dataf[DataFrames.isna(dataf[i]), i] = "NONE"
        elseif typeof(dataf[i]) == DataArrays.DataArray{UTF8String,1}
              dataf[DataFrames.isna(dataf[i]), i] = "NONE"
      end
        if sum(size(dataf[DataFrames.isna(dataf[i]), i]))>0
        print(i, "\n")
      end
    end
    data = convert(Matrix, dataf);
    dataf = 0; #set to zero to clear out.

      medicaidmodcoeffs = DataFrames.readtable(pathpeople*"TX 2005 Medicaid Model.csv", header = true);
      global const medicaiddistance_c = medicaidmodcoeffs[1, 2]
      global const medicaiddistsq_c = medicaidmodcoeffs[2, 2]
      global const medicaidneoint_c = medicaidmodcoeffs[3, 2]
      global const medicaidsoloint_c = medicaidmodcoeffs[4, 2]
      global const medicaidclosest_c = medicaidmodcoeffs[5, 2]
      global const medicaiddistbed_c = medicaidmodcoeffs[6, 2]
      medicaiddemandmodelparameters = [medicaiddistance_c medicaiddistsq_c medicaidneoint_c medicaidsoloint_c medicaidclosest_c medicaiddistbed_c]

      privatemodcoeffs = DataFrames.readtable(pathpeople*"TX 2005 Private Ins Model.csv", header = true);
      global const privatedistance_c = privatemodcoeffs[1, 2]
      global const privatedistsq_c = privatemodcoeffs[2, 2]
      global const privateneoint_c = privatemodcoeffs[3, 2]
      global const privatesoloint_c = privatemodcoeffs[4, 2]
      global const privateclosest_c = privatemodcoeffs[5, 2]
      global const privatedistbed_c = privatemodcoeffs[6, 2]
      privatedemandmodelparameters = [privatedistance_c privatedistsq_c privateneoint_c privatesoloint_c privateclosest_c privatedistbed_c]

    println("Importing Medicaid Patients") # use the infants only.
    medicaid = DataFrames.readtable(pathpeople*"TX 2005 Medicaid Individual Choices.csv", header = true);
    #people = DataFrames.readtable(pathpeople*"TX 2005 1 Individual Choices.csv", header = true); #smaller version for testing.
    for i in names(medicaid)
      if ( typeof(medicaid[i]) == DataArrays.DataArray{Float64,1} )
        medicaid[DataFrames.isna(medicaid[i]), i] = 0
      elseif (typeof(medicaid[i]) == DataArrays.DataArray{Int64,1})
        medicaid[DataFrames.isna(medicaid[i]), i] = 0
      elseif typeof(medicaid[i]) == DataArrays.DataArray{ByteString,1}
        # A dumb way to make sure no one chooses a missing facility: set covariate values to large numbers
        # with opposite signs of the corresponding coefficients from modelparameters.
        # This does that by looking at missing NAMES, not fids.
        medicaid[DataFrames.isna(medicaid[i]), medicaid.colindex.lookup[i]+2] = -sign(medicaidneoint_c)*99
        medicaid[DataFrames.isna(medicaid[i]), medicaid.colindex.lookup[i]+8] = -sign(medicaidsoloint_c)*99
        medicaid[DataFrames.isna(medicaid[i]), i] = "NONE"
      elseif typeof(medicaid[i]) == DataArrays.DataArray{UTF8String,1}
        medicaid[DataFrames.isna(medicaid[i]), medicaid.colindex.lookup[i]+2] = -sign(medicaidneoint_c)*99
        medicaid[DataFrames.isna(medicaid[i]), medicaid.colindex.lookup[i]+8] = -sign(medicaidsoloint_c)*99
        medicaid[DataFrames.isna(medicaid[i]), i] = "NONE"
      end
      if sum(size(medicaid[DataFrames.isna(medicaid[i]), i]))>0
        println(i)
      end
    end
    pmedicaid = convert(Matrix, medicaid);
    medicaid = 0; # DataFrame not used - set to 0 and clear out.
    for i =1:size(pmedicaid, 2)
      if (typeof(pmedicaid[2,i])==UTF8String) | (typeof(pmedicaid[2,i])==ASCIIString)
        #      print(i, "\n")
        pmedicaid[:,i] = "0"
        pmedicaid[:,i] = map(x->parse(Float64, x), pmedicaid[:,i])
      end
    end
    # Note this change - I don't think there's anything that requires 64 bits.
    pmedicaid = convert(Array{Float32, 2}, pmedicaid)
    println("Size of Medicaid, ", size(pmedicaid))

    println("Importing Privately Insured Patients") #use the infants only.
    pinsure = DataFrames.readtable(pathpeople*"TX 2005 Private Ins Individual Choices.csv", header = true);
    for i in names(pinsure)
      if ( typeof(pinsure[i]) == DataArrays.DataArray{Float64,1} )
        pinsure[DataFrames.isna(pinsure[i]), i] = 0
      elseif (typeof(pinsure[i]) == DataArrays.DataArray{Int64,1})
        pinsure[DataFrames.isna(pinsure[i]), i] = 0
      elseif typeof(pinsure[i]) == DataArrays.DataArray{ByteString,1}
        # A dumb way to make sure no one chooses a missing facility: set covariate values to large numbers
        # with opposite signs of the corresponding coefficients from modelparameters.
        # This does that by looking at missing NAMES, not fids.
        pinsure[DataFrames.isna(pinsure[i]), pinsure.colindex.lookup[i]+2] = -sign(privateneoint_c)*99
        pinsure[DataFrames.isna(pinsure[i]), pinsure.colindex.lookup[i]+8] = -sign(privatesoloint_c)*99
        pinsure[DataFrames.isna(pinsure[i]), i] = "NONE"
      elseif typeof(pinsure[i]) == DataArrays.DataArray{UTF8String,1}
        pinsure[DataFrames.isna(pinsure[i]), pinsure.colindex.lookup[i]+2] = -sign(privateneoint_c)*99
        pinsure[DataFrames.isna(pinsure[i]), pinsure.colindex.lookup[i]+8] = -sign(privatesoloint_c)*99
        pinsure[DataFrames.isna(pinsure[i]), i] = "NONE"
      end
      if sum(size(pinsure[DataFrames.isna(pinsure[i]), i]))>0
        println(i)
      end
    end
     pinsured = convert(Matrix,pinsure);
     pinsure= 0; # DataFrame not used - set to 0 and clear out.
    for i =1:size(pinsured, 2)
      if (typeof(pinsured[2,i])==UTF8String) | (typeof(pinsured[2,i])==ASCIIString)
  #      print(i, "\n")
        pinsured[:,i] = "0"
        pinsured[:,i] = map(x->parse(Float64, x), pinsured[:,i])
      end
    end
    # Note this change - I don't think there's anything that requires 64 bits.
     pinsured= convert(Array{Float32, 2}, pinsured)
     println("Size of Privately Insured, ", size(pinsured))
     # TX zip codes:
     TXzps = DataFrames.readtable(pathprograms*"TXzipsonly.csv", header = false)
     TXzips = convert(Vector, TXzps[:x1])
     # DRG codes:
     # These are for infants only.
     DRGs = [385 386 387 388 389 390 391]

    global const fipscodeloc = 78; # this is for hospital data, here as "data"
    global const yearloc = 75; # this also for hospital data, here imported as "data"
    global const fidloc = 74; # Also for hospital data, here as "data"
    global const idloc = 1; # Also for Hospital data, here as "data"
    # Collect FIDs
    txfd = DataFrames.readtable(pathprograms*"TXfidsonly.csv", header = true)
    allfids = convert(Vector, txfd[:fid])

    # WTP values -
    txwp = DataFrames.readtable(pathprograms*"WTPTemplate.csv", header = true)
    WTP = convert(Matrix, txwp)
  #  global const nsims = 500;
  #  global const dirs = pwd() # present working directory path
end # of "begin" block


#
# Test functions:
println("Locations - are they correct?")
FindCorrect(pinsured)
println("Testing Demand")
test = DetUtil(pinsured, privatedemandmodelparameters)
println("Testing WTP Computation")
testWTP = ComputeWTP(test)
println("WTP Mapping")
mappedWTP = MapWTP(testWTP)
println("Check Return WTP")
sumWTP = ReturnWTP(mappedWTP)
print("Testing Simulator \n")
a1 = Simulator(data, pinsured, privatedemandmodelparameters, pmedicaid, medicaiddemandmodelparameters, 2003, 48027 ; T = 1) #tests logitest, DemandModel, distance
print("Testing PerturbSimulator \n")
 b1 = PerturbSimulator(data, pinsured, privatedemandmodelparameters, pmedicaid, medicaiddemandmodelparameters, 2003, 48027, 273410; T = 1)
print("Testing logitest \n")
logitest( (1, 0), 2, 3, 1, zeros(9))
print("Testing DemandModel \n")
DemandModel(pinsured, privatedemandmodelparameters, [99999 1 0 120 32.96  -96.8385], true) #implicitly tests EntrantsU as well
print("Testing DynamicValue \n")
DynamicValue(b1, b1[1,1])

#Test parallel
print("Testing Remote Call \n")
#p1 = remotecall_fetch(lis[1], DemandModel, pinsured, privatedemandmodelparameters, Array{Float64,2}())
#p2 = remotecall_fetch(lis[2], DemandModel, pmedicaid, medicaiddemandmodelparameters, Array{Float64, 2}())

#include(pathprograms*"Main.jl")





###
