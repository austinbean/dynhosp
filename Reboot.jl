


#begin
  ##  Figure out which machine I'm on

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
    global pathdata = "/home1/04179/abean/dynhosp/"
    global pathpeople = "/home1/04179/abean/dynhosp/"
    global pathprograms = "/home1/04179/abean/dynhosp/"
  end
  ###### Import the hospital data and convert to a matrix -


    # TODO -

  println("Importing Hosp Data")
    dataf = DataFrames.readtable(pathdata*"TX Transition Probabilities.csv", header = true);
    for i in names(dataf)
        if ( typeof(dataf[i]) == DataArrays.DataArray{Float64,1} )
            dataf[DataFrames.isna(dataf[i]), i] = 0
        elseif (typeof(dataf[i]) == DataArrays.DataArray{Int64,1})
            dataf[DataFrames.isna(dataf[i]), i] = 0
        elseif typeof(dataf[i]) == DataArrays.DataArray{String,1} #Changed for 0.5
            dataf[DataFrames.isna(dataf[i]), i] = "NONE"
        elseif typeof(dataf[i]) == DataArrays.DataArray{String,1}
              dataf[DataFrames.isna(dataf[i]), i] = "NONE"
      end
        if sum(size(dataf[DataFrames.isna(dataf[i]), i]))>0
        print(i, "\n")
      end
    end
    data = convert(Matrix, dataf);
    dataf = 0; #set to zero to clear out.


#DONE 
      medcoeffs, medcoeffsnames = readcsv(pathpeople*"TX 2005 Medicaid Model.csv", header = true);

      global const medicaiddistance_c = medcoeffs[1,2] 
      global const medicaiddistsq_c = medcoeffs[2,2]
      global const medicaidneoint_c = medcoeffs[3,2]
      global const medicaidsoloint_c = medcoeffs[4,2]
      global const medicaidclosest_c = medcoeffs[5,2]
      global const medicaiddistbed_c = medcoeffs[6,2]
      medicaiddemandmodelparameters = [medicaiddistance_c medicaiddistsq_c medicaidneoint_c medicaidsoloint_c medicaidclosest_c medicaiddistbed_c]

      privcoeffs, privcoeffsnames = readcsv(pathpeople*"TX 2005 Private Ins Model.csv", header = true)

      global const privatedistance_c = privcoeffs[1,2]
      global const privatedistsq_c = privcoeffs[2,2]
      global const privateneoint_c = privcoeffs[3,2]
      global const privatesoloint_c = privcoeffs[4,2]
      global const privateclosest_c = privcoeffs[5,2]
      global const privatedistbed_c = privcoeffs[6,2]
      privatedemandmodelparameters = [privatedistance_c privatedistsq_c privateneoint_c privatesoloint_c privateclosest_c privatedistbed_c]



#DONE 
    println("Importing Medicaid Patients") # use the infants only.
    pmedicaid, medicaidnames = readcsv(pathpeople*"TX 2005 Medicaid Individual Choices.csv", header = true)

    for i = 1:size(pinsured,1)
        for j = 1:size(pinsured,2)
            if typeof(pinsured[i,j])<:AbstractString
                if (pinsured[i,j]=="")
                    pinsured[i,j] = 0 # = 0
                end 
            end
        end 
    end 

    pmedicaid = convert(Array{Float32, 2}, pmedicaid)


#DONE 
    println("Importing Privately Insured Patients") #use the infants only.
    
    pinsured, pinsurednames = readcsv(pathpeople*"TX 2005 Private Ins Individual Choices.csv", header = true)


    for i = 1:size(pinsured,1)
        for j = 1:size(pinsured,2)
            if typeof(pinsured[i,j])<:AbstractString
                if (pinsured[i,j]=="")
                    pinsured[i,j] = 0 # = 0
                end 
            end
        end 
    end 
    # Note this change - I don't think there's anything that requires 64 bits.
     pinsured= convert(Array{Float32, 2}, pinsured)


  # DONE 

     zp = readcsv(pathprograms*"TXzipsonly.csv", header = false)
     TXzips = convert(Array{Int64,1}, zp[:,1])



     # DRG codes:
     # These are for infants only.
     DRGs = [385 386 387 388 389 390 391]

    global const fipscodeloc = 78; # this is for hospital data, here as "data"
    global const yearloc = 75; # this also for hospital data, here imported as "data"
    global const fidloc = 74; # Also for hospital data, here as "data"
    global const idloc = 1; # Also for Hospital data, here as "data"
    # Collect FIDs
    #DONE 
    allfids, txfnames = readcsv(pathprograms*"TXfidsonly.csv", header = true)
    allfids = convert(Vector, allfids)



#DONE
  zips = readcsv(pathprograms*"TXzipsonly.csv", header = false);
  zips = convert(Array, zp[:,1]);




  # the choices file above is not correct.  Load TX All Distances.csv instead:
  # Get this from hospitaldistancepair.py and TX Hospital Sets.do and hosplatlong.py
  #TODO - replace all NA's with zeros.  Change all calls to ProjectModule.alldists which are dataframes to array.
  # 02/25/2017 - I am still loading both of these files even though they are redundant and it seems the second is wrong.


  # TODO
  alldists = DataFrames.readtable(pathdata*"TX Zip All Hospital Distances.csv", header = true);
  for n in alldists.colindex.names
    if typeof(alldists[n]) == DataArrays.DataArray{String,1}
      alldists[ isna(alldists[n]),n] = "missing"
      println("changed") # there are 3 places only where this is changed.  
    else
      alldists[ isna(alldists[n]), n] = 0
    end
  end
  alldists = convert(Array, alldists)

  alld, alldlabs = readcsv(pathdata*"TX Zip All Hospital Distances.csv", header = true);

  for i = 1:size(alld,1)
    for j = 1:size(alld,2)
      if alld[i,j] == ""
        alld[i,j] = 0
      end 
    end 
  end 


diffs = Set()
difflocs = Set()
diffsandlocs = Set()

for i = 1:size(alld,1)
    for j = 1:size(alld,2)
        if alld[i,j]!=alldists[i,j]
            push!(diffs, (alld[i,j], alldists[i,j]))
            push!(difflocs, (i,j))
            push!(diffsandlocs, (i,j, alld[i,j], alldists[i,j]))
        end 
    end 
end 










  choices = DataFrames.readtable(pathdata*"TX Zip Code Choice Sets.csv", header = true);
  for el in choices.colindex.names
    #println(typeof(choices[el]), "  ", el)
    if (typeof(choices[el]) == DataArrays.DataArray{Int64,1})|(typeof(choices[el]) == DataArrays.DataArray{Float64,1})
      choices[isna(choices[:,el]) , el] = 0
    elseif (typeof(choices[el]) == DataArrays.DataArray{String,1})
      choices[isna(choices[:,el]), el] = "Missing"
    end
  end

  choices = convert(Array{Any, 2}, choices);

  # TODO  - Data
  fips = convert(Array{Int64}, union(unique(data[:,78]), unique(alldists[:,7])))
  data05 = 0 # data[(data[:,75].==2005), :] ; # NB - killing this because it has bad data and I want to find what is using it and change it.

  #DONE  - 

  bwp, bwplabel = readcsv(pathdata*"2005 Birth Weight Probabilities.csv", header = true);

  np, nplabel = readcsv(pathdata*"2005 NICU Admission Probabilities.csv", header = true);

  weightprobs = zeros(size(bwprobs,1), 2)
  nicuprobs = zeros(size(naprobs,1), 2)
  for el = 1:size(bwp,1)
    weightprobs[el, 1] = el
    nicuprobs[el,1] = el
    weightprobs[el, 2] = bwp[el,2]
    nicuprobs[el,2] = np[el,2]
  end

  bwp = 0;
  np = 0;


#end # of "begin" block





###
