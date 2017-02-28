


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
  elseif (dir=="/work/04179/abean/dynhosp")
    global pathdata = "/work/04179/abean/dynhosp/"
    global pathpeople = "/work/04179/abean/dynhosp/"
    global pathprogram = "/work/04179/abean/dynhosp/"
  end
  ###### Import the hospital data and convert to a matrix -


# Data Concerning Hospital Choices.

  println("Importing Hosp Data")

  data, datnames = readcsv(pathdata*"TX Transition Probabilities.csv", header = true, comments = false); # the use of the \# sign made the csvreader think there were comments.  
  for j = 1:size(data,2) # covers the first row separately
    if (data[1,j] == "")
      if (typeof(data[1+1,j])<:AbstractString)&(length(data[1+1,j])>0)
          data[1, j] = "NONE" # this branch never triggers for the array dat.
        else 
          data[1,j] = 0
        end 
    end 
  end
  for i = 2:size(data,1)
    for j = 1:size(data,2)
      if (data[i,j] == "")
        if (typeof(data[i-1,j])<:AbstractString)&(length(data[i-1,j])>0)
          if (data[i-1,j]=="0-5 Miles")||(data[i-1,j]=="5-15 Miles")||(data[i-1,j]=="15-25 Miles")
            data[i, j] = data[i-1,j] # matches previous line.
          else 
            data[i,j]="NONE"
          end 
        else 
          data[i,j] = 0
        end 
      end 
    end 
  end

#Model Coefficients for Privately Insured and Medicaid Patients.  
      medcoeffs, medcoeffsnames = readcsv(pathpeople*"TX 2005 Medicaid Model.csv", header = true);

      global const medicaiddistance_c = medcoeffs[1,2] 
      global const medicaiddistsq_c = medcoeffs[2,2]
      global const medicaidneoint_c = medcoeffs[3,2]
      global const medicaidsoloint_c = medcoeffs[4,2]
      global const medicaidclosest_c = medcoeffs[5,2]
      global const medicaiddistbed_c = medcoeffs[6,2]
      medicaiddemandmodelparameters = [medicaiddistance_c medicaiddistsq_c medicaidneoint_c medicaidsoloint_c medicaidclosest_c medicaiddistbed_c]

      privcoeffs, privcoeffsnames = readcsv(pathpeople*"TX 2005 Private Ins Model.csv", header = true);

      global const privatedistance_c = privcoeffs[1,2]
      global const privatedistsq_c = privcoeffs[2,2]
      global const privateneoint_c = privcoeffs[3,2]
      global const privatesoloint_c = privcoeffs[4,2]
      global const privateclosest_c = privcoeffs[5,2]
      global const privatedistbed_c = privcoeffs[6,2]
      privatedemandmodelparameters = [privatedistance_c privatedistsq_c privateneoint_c privatesoloint_c privateclosest_c privatedistbed_c]

#Medicaid Patients. 
    println("Importing Medicaid Patients") # use the infants only.
    pmedicaid, medicaidnames = readcsv(pathpeople*"TX 2005 Medicaid Individual Choices.csv", header = true);

    for i = 1:size(pmedicaid,1)
        for j = 1:size(pmedicaid,2)
            if typeof(pmedicaid[i,j])<:AbstractString
                if (pmedicaid[i,j]=="")
                    pmedicaid[i,j] = 0 # = 0
                end 
            end
        end 
    end 

    pmedicaid = convert(Array{Float32, 2}, pmedicaid);

#Privately Insured Patients.  
    println("Importing Privately Insured Patients") #use the infants only.
    
    pinsured, pinsurednames = readcsv(pathpeople*"TX 2005 Private Ins Individual Choices.csv", header = true);


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
     pinsured= convert(Array{Float32, 2}, pinsured);


  # Zip Codes - a list of all of them.   
     zp = readcsv(pathprograms*"TXzipsonly.csv", header = false)
     TXzips = convert(Array{Int64,1}, zp[:,1])
  #Zip Codes - a list of all of them 
    zips = readcsv(pathprograms*"TXzipsonly.csv", header = false);
    zips = convert(Array, zp[:,1]);



     # DRG codes:
     # These are for infants only.
     DRGs = [385 386 387 388 389 390 391]

    global const fipscodeloc = 78; # this is for hospital data, here as "data"
    global const yearloc = 75; # this also for hospital data, here imported as "data"
    global const fidloc = 74; # Also for hospital data, here as "data"
    global const idloc = 1; # Also for Hospital data, here as "data"

    # Collect FIDs
    allfids, txfnames = readcsv(pathprograms*"TXfidsonly.csv", header = true)
    allfids = vec(convert(Array{Int64,2}, allfids))

  # Get this from hospitaldistancepair.py and TX Hospital Sets.do and hosplatlong.py
  # 02/25/2017 - I am still loading both of these files even though they are redundant and it seems the second is wrong.


# Distances 
  alldists, alldlabs = readcsv(pathdata*"TX Zip All Hospital Distances.csv", header = true);

  for j = 1:size(alldists,2) # covers the first row separately
    if (alldists[1,j] == "")
      alldists[1,j] = 0.0
    end 
  end
  for i = 2:size(alldists,1)
    for j = 1:size(alldists,2)
      if (alldists[i,j] == "")
        if (typeof(alldists[i-1,j])<:AbstractString)&(length(alldists[i-1,j])>0)
          alldists[i, j] = "missing" # matches the previous importation strategy.
        else 
          alldists[i,j] = 0
        end 
      end 
    end 
  end 

#Zip Code Choice Sets 
  choices, chnames = readcsv(pathdata*"TX Zip Code Choice Sets.csv", header = true);

  for j = 1:size(choices,2) # covers the first row separately
    if (choices[1,j] == "")
      if (typeof(choices[1+1,j])<:AbstractString)&(length(choices[1+1,j])>0)
          choices[1, j] = "Missing" # matches the previous importation strategy.
        else 
          choices[1,j] = 0
        end 
    end 
  end
  for i = 2:size(choices,1)
    for j = 1:size(choices,2)
      if (choices[i,j] == "")
        if (typeof(choices[i-1,j])<:AbstractString)&(length(choices[i-1,j])>0)
          choices[i, j] = "Missing" # matches the previous importation strategy.
        else 
          choices[i,j] = 0
        end 
      end 
    end 
  end 



  fips = convert(Array{Int64}, union(unique(data[:,78]), unique(alldists[:,7])))
  data05 = 0 # data[(data[:,75].==2005), :] ; # NB - killing this because it has bad data and I want to find what is using it and change it.

# Birthweight and NICU Admission Probs 
  bwp, bwplabel = readcsv(pathdata*"2005 Birth Weight Probabilities.csv", header = true);

  np, nplabel = readcsv(pathdata*"2005 NICU Admission Probabilities.csv", header = true);

  weightprobs = zeros(size(bwp,1), 2)
  nicuprobs = zeros(size(np,1), 2)
  for el = 1:size(bwp,1)
    weightprobs[el, 1] = el
    nicuprobs[el,1] = el
    weightprobs[el, 2] = bwp[el,2]
    nicuprobs[el,2] = np[el,2]
  end

  bwp = 0;
  np = 0;


#end # of "begin" block


#=
# checking differences between old and new import strategies.
diffs = Set()
difflocs = Set()
diffsandlocs = Set()

for i = 1:size(dat,1)
    for j = 1:size(dat,2)
        if dat[i,j]!=data[i,j]
            push!(diffs, (dat[i,j], data[i,j]))
            push!(difflocs, (i,j))
            push!(diffsandlocs, (i,j, dat[i,j], data[i,j]))
        end 
    end 
end 
=#



###
