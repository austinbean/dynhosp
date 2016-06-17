# Try to do this Parallelizing with main function.

# Parallelize over nsims, I think.

lis = addprocs()

push!(LOAD_PATH, "/Users/austinbean/Desktop/dynhosp")
push!(LOAD_PATH, "/dynhosp/dynhosp")
push!(LOAD_PATH, "/home/ubuntu/dynhosp/")







 using ProjectModule
@everywhere using DataFrames
@everywhere using Distributions

@everywhere begin
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
  end
  # Import the hospital data and convert to a matrix -
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


      # Import the model coefficients
      modcoeffs = DataFrames.readtable(pathpeople*"TX 2005 Model.csv", header = true);

      global const distance_c = modcoeffs[1, 2]
      global const distsq_c = modcoeffs[2, 2]
      global const neoint_c = modcoeffs[3, 2]
      global const soloint_c = modcoeffs[4, 2]
      global const closest_c = modcoeffs[5, 2]
      global const distbed_c = modcoeffs[6, 2]

      demandmodelparameters = [distance_c distsq_c neoint_c soloint_c closest_c distbed_c]



  # Import the people and convert that data to a matrix
  #  people = DataFrames.readtable(pathpeople*"TX 2005 Individual Choices.csv", header = true);
    people = DataFrames.readtable(pathpeople*"TX 2005 1 Individual Choices.csv", header = true); #smaller version for testing.

        for i in names(people)
          if ( typeof(people[i]) == DataArrays.DataArray{Float64,1} )
            people[DataFrames.isna(people[i]), i] = 0
          elseif (typeof(people[i]) == DataArrays.DataArray{Int64,1})
            people[DataFrames.isna(people[i]), i] = 0
          elseif typeof(people[i]) == DataArrays.DataArray{ByteString,1}
            # A dumb way to make sure no one chooses a missing facility: set covariate values to large numbers
            # with opposite signs of the corresponding coefficients from modelparameters.
            # This does that by looking at missing NAMES, not fids.
            people[DataFrames.isna(people[i]), people.colindex.lookup[i]+2] = -sign(neoint_c)*99
            people[DataFrames.isna(people[i]), people.colindex.lookup[i]+8] = -sign(soloint_c)*99
            people[DataFrames.isna(people[i]), i] = "NONE"
          elseif typeof(people[i]) == DataArrays.DataArray{UTF8String,1}
            people[DataFrames.isna(people[i]), people.colindex.lookup[i]+2] = -sign(neoint_c)*99
            people[DataFrames.isna(people[i]), people.colindex.lookup[i]+8] = -sign(soloint_c)*99
            people[DataFrames.isna(people[i]), i] = "NONE"
          end
          if sum(size(people[DataFrames.isna(people[i]), i]))>0
            print(i, "\n")
          end
        end

        peoples = convert(Matrix, people);
        people = 0; # DataFrame not used - set to 0 and clear out.


        for i =1:size(peoples, 2)
          if (typeof(peoples[2,i])==UTF8String) | (typeof(peoples[2,i])==ASCIIString)
            print(i, "\n")
            peoples[:,i] = "0"
            peoples[:,i] = map(x->parse(Float64, x), peoples[:,i])
          end
        end
        peoples = convert(Array{Float64, 2}, peoples)

        global const fipscodeloc = 78;
        global const yearloc = 75;
        global const fidloc = 74;
        global const idloc = 1;

      #  global const dirs = pwd() # present working directory path


  end # of "begin" block

    @everywhere begin
    yearins = [ [x; findfirst(data[:,fipscodeloc], x); findlast(data[:,fipscodeloc], x ); unique( data[findfirst(data[:,fipscodeloc], x):findlast(data[:,fipscodeloc], x ) , yearloc]  ) ] for x in unique(data[:,fipscodeloc])  ]
    mkt_fips = yearins[11][1]
    year = 2005
    fids = convert(Array{Int64, 1}, sort!(unique(data[(data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year),fidloc])))
    end


  # Parallel Function Tests:

  p1 = remotecall_fetch(lis[1], DemandModel, peoples, demandmodelparameters, Array{Float64,2}())
  p2 = remotecall_fetch(lis[2], Simulator, data, peoples, year, mkt_fips, demandmodelparameters)


@everywhere function ParMainfun(dataf::Matrix, people::Matrix, mkt_fips::Int64, year::Int64, modelparameters::Array{Float64, 2}, fids::Array{Int64};  idloc = 1, nsims = 3, npers = 10)
			 # returns a dataframe unless converted
      numfids = maximum(size(fids))
      outp = zeros(numfids, 183)
      for i = 1:numfids
        outp[i,1] = mkt_fips
        outp[i,2] = fids[i]
        outp[i,3] = year
      end
      states =
			dataf = dataf[(dataf[:,idloc].>= 0), :];
      states = ProjectModule.Simulator(dataf, people, year, mkt_fips, modelparameters; T = npers)
  		for f in 1:numfids
          pfid = fids[f]
          pfid_f = convert(Float64, pfid)
          # Note that it is important that all user-defined functions in this function say ProjectModule.(fn) - these are imported by the module but not brought into scope.
          # This probably won't be necessary if this is a function in its own file and imported directly via the module.  (The other two functions would surely cause problems if that were necessary.)
          neq_change, neq_val = ProjectModule.DynamicValue(ProjectModule.PerturbSimulator(dataf, people, year, mkt_fips, modelparameters, pfid; disturb = 0.01, T = npers), pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
          eq_change, eq_val = ProjectModule.DynamicValue(states, pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
          # Abandon Entrants again.
          dataf = dataf[(dataf[:,idloc].>= 0), :];
          outp[f,4:end] += [eq_val eq_change neq_val neq_change]
      end
  return outp
end

p2 = remotecall_fetch(lis[2], ParMainfun, data, peoples, mkt_fips, year, demandmodelparameters, [fids[1]])


# The parallelized version - this works.
# Extremely easy conditions: 5.187119 seconds (18.76 M allocations: 1.135 GB, 2.13% gc time)
# Slightly more demand: nsims = 3, npers = 10 - 13.411527 seconds (18.77 M allocations: 1.135 GB, 1.22% gc time)
@parallel (+) for i = 1:3
  ParMainfun(data, peoples, mkt_fips, year, demandmodelparameters, [fids[1]])
end

# For comparison:
# nsims = 3, npers = 10 - 24.139081 seconds (65.39 M allocations: 8.865 GB, 17.99% gc time)
a2 = Mainfun(data, peoples, mkt_fips, year, demandmodelparameters, [fids[1]]; nsims = 4, npers = 10)


@parallel (+) for i=1:4
  ParMainfun(data, peoples, mkt_fips, year, demandmodelparameters, fids;)
end



print(outp)
