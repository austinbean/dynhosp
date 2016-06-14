# Main

#=
Version Information:
Julia Version 0.4.3
Commit a2f713d (2016-01-12 21:37 UTC)
Platform Info:
  System: Darwin (x86_64-apple-darwin13.4.0)
  CPU: Intel(R) Core(TM) i7-5557U CPU @ 3.10GHz
  WORD_SIZE: 64
  BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
  LAPACK: libopenblas64_
  LIBM: libopenlibm
  LLVM: libLLVM-3.3
include("/Users/austinbean/Desktop/dynhosp/Main.jl")
=#

# Hypothetical:


push!(LOAD_PATH, "/Users/austinbean/Desktop/dynhosp")
push!(LOAD_PATH, "/dynhosp/dynhosp")
push!(LOAD_PATH, "/home/ubuntu/dynhosp/")

using Reboot.jl



### Collect Basic Information ###

# locates starting and ending points of all separate fipscode indices

allindices = [ (x, findfirst(data[:,fipscodeloc], x), findlast(data[:,fipscodeloc], x)) for x in unique(data[:,fipscodeloc]) ];


# Next one stores all the years appearing in each fips code
yearins = [ [x; findfirst(data[:,fipscodeloc], x); findlast(data[:,fipscodeloc], x ); unique( data[findfirst(data[:,fipscodeloc], x):findlast(data[:,fipscodeloc], x ) , yearloc]  ) ] for x in unique(data[:,fipscodeloc])  ];



# Action codes:
#=
1 "To 3 from 1"
2 "To 2 from 1"
3 "To 2 from 3"
4 "To 1 from 3"
5 "To 1 from 2"
6 "To 3 from 2"
7 "Enter at 1"
8 "Enter at 2"
9 "Enter at 3"
10 "Do Nothing"
11 "Exit"
=#



#=
Future direction - write this to start where it left off.
# Start with some smaller markets first:

=#
monopoly = Array{Int64}(0)
duopoly = Array{Int64}(0)
triopoly = Array{Int64}(0)
tetrapoly = Array{Int64}(0)
nopoly = Array{Int64}(0)

for el in yearins
  unqfids = [x for x in unique(data[el[2]:el[3],fidloc])]
  if size(unqfids)[1] == 1
  #  print("Fipscode Monopoly: ", unique(dataf[el[2]:el[3], :fipscode]), "\n")
    push!(monopoly, data[el[3], fipscodeloc])
  elseif size(unqfids)[1] == 2
  #  print("Fipscode Duopoly: ", unique(dataf[el[2]:el[3], :fipscode]), "\n")
    push!(duopoly, data[el[3], fipscodeloc])
  elseif size(unqfids)[1] == 3
  #  print("Fipscode Triopoly: ", unique(dataf[el[2]:el[3], :fipscode]), "\n")
    push!(triopoly, data[el[3], fipscodeloc])
  elseif size(unqfids)[1] == 4
  #  print("Fipscode Tetrapoly: ", unique(dataf[el[2]:el[3], :fipscode]), "\n")
    push!(tetrapoly, data[el[3], fipscodeloc])
  elseif size(unqfids)[1] > 4
  #  print("Fipscode N-opoly: ", unique(dataf[el[2]:el[3], :fipscode]), "\n")
  #  print("Fipscode Hospitals: ", size(unqfids)[1], "\n")
    push!(nopoly, data[el[3], fipscodeloc])
  end
end


container = zeros(1, 183)

# Open the existing saved data:
fout1 = DataFrames.readtable(pathprograms*"simulationresults.csv")
donefips  = [x for x in unique(fout1[!(DataFrames.isna(fout1[:fipscode])),:fipscode])] # take codes NOT done before
#donefips = [] # set this to be empty until the simulation gets working.

entryprobs = [0.9895, 0.008, 0.0005, 0.002]


colnames = Array{Symbol}(:0)
push!(colnames, :fipscode)
push!(colnames, :fid)
push!(colnames, :year)
for elem in ["EQ", "NEQ"]
  for j = 1:3
    for k = 0:25
      push!(colnames, parse("$elem"*"Lev$j"*"Comp$k"))
    end
  end
  for x in [1 2 3]
    for y in [1 2 3 "EX"]
      if x != y
        push!(colnames, parse("$elem"*"Trans$x$y"))
      end
    end
  end
  push!(colnames, parse("$elem"*"Enter1"))
  push!(colnames, parse("$elem"*"Enter2"))
  push!(colnames, parse("$elem"*"Enter3"))
end

timestamps = Array{Any}(0)

#=
for y in 1#:size(duopoly,1)
    mkt_fips = duopoly[y][1]
    crtime = now()
    timestr = Dates.format(crtime, "yyyy-mm-dd HH:MM:ss")
    push!(timestamps, (mkt_fips, "begin", timestr))
    if !(in(mkt_fips, donefips)) # this is going to do new fipscodes only
      print("Market FIPS Code ", mkt_fips, "\n")
      	for year in [ 2005 ]   #yearins[y][4:end] # can do all years or several.
          #dat = deepcopy(data);
          #print("size of dat", size(dat), "\n")
          fids =  sort!(convert(Array{Int64}, unique(data[(data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year),fidloc])))
          numfids = size(fids,1)
          container = [container; Mainfun(data, peoples, mkt_fips, year, demandmodelparameters, fids; nsims = 52, npers = 2)]
      end #Note - rewrite first line of state history back to peoples.
    #  dat = 0
    end
    # Record what time certain ends happened.
    outf = open(pathprograms*"timer.txt", "w")
    crtime = now()
    timestr = Dates.format(crtime, "yyyy-mm-dd HH:MM:ss")
    push!(timestamps, (mkt_fips, "end", timestr))
    #print(timestamps, "\n")
    writedlm(outf, timestamps)
    close(outf)
    # Write out temporary files in case the process fails later.
    tout = convert(DataFrame, container);
    names!(tout, colnames)
    tout = tout[ tout[:fipscode].>0, :]
    DataFrames.writetable(pathprograms*"temp_results_$mkt_fips.csv", tout) # no slash.
end

=#

for y in 1#:size(duopoly,1)
    mkt_fips = duopoly[y][1]
    crtime = now()
    timestr = Dates.format(crtime, "yyyy-mm-dd HH:MM:ss")
    push!(timestamps, (mkt_fips, "begin", timestr))
    if !(in(mkt_fips, donefips)) # this is going to do new fipscodes only
      print("Market FIPS Code ", mkt_fips, "\n")
      	for year in [ 2005 ]   #yearins[y][4:end] # can do all years or several.
          #dat = deepcopy(data);
          #print("size of dat", size(dat), "\n")
          fids =  sort!(convert(Array{Int64}, unique(data[(data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year),fidloc])))
          numfids = maximum(size(fids))
          outp = zeros(numfids, 183)
          for i = 1:numfids
            outp[i,1] = mkt_fips
            outp[i,2] = fids[i]
            outp[i,3] = year
          end
          for j in 1:51  # NUMBER OF SIMULATIONS
          #  print("Sim Round, " , j, "\n")
      			#Arguments: Simulator(data::Matrix, peoplesub::Matrix, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2};  T = 100, sim_start = 2, fields = 7, neighbors_start = 108, entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002])
            if j%50 == 0
              print("Simulation Round, ", j, "  mkt: ", mkt_fips, " ", year, " ", "\n")
            end
            states = Simulator(data, peoples, year, mkt_fips, demandmodelparameters; T = 2) # T is PERIODS SIMULATED
            # Non-equilibrium Play -
      			# Entrants in dataframe now tagged with negative ID's.  Remake to remove them:
      			data = data[(data[:,idloc].>= 0), :];
        		for f in 1:numfids
                pfid = fids[f]
            #    print("Perturbing Fid: ", pfid, "\n")
          			#Arguments: function PerturbSimulator(data::Matrix, peoplesub::Matrix, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2}, pfid::Int64; entryprobs = [0.9895, 0.008, 0.0005, 0.002], entrants = [0, 1, 2, 3], disturb = 0.05,  T = 100, sim_start = 2, fields = 7, neighbors_start = 108)
                perturbed_history = PerturbSimulator(data, peoples, year, mkt_fips, demandmodelparameters, pfid; disturb = 0.01, T = 2)  # T is PERIODS SIMULATED

          			# Here apply DynamicValue to the result of the simulations
          			# DynamicValue(state_history::Array, fac_fid::Float64; α₂ = 0.07, α₃ = 0.13, pat_types = 1, β = 0.95, max_hosp = 25)
          			# output is in format: facility changes record, per-period visits record.
          			pfid_f = convert(Float64, pfid)
          			eq_change, eq_val  = DynamicValue(states, pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
          			neq_change, neq_val = DynamicValue(perturbed_history, pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
                outp[f,4:end] += [eq_val eq_change neq_val neq_change]
                # Abandon Entrants again.
                data = data[(data[:,idloc].>= 0), :];
            end
        end
        container = [container; outp]
      end #Note - rewrite first line of state history back to peoples.
    #  dat = 0
    end
    # Record what time certain ends happened.
    outf = open(pathprograms*"timer.txt", "w")
    crtime = now()
    timestr = Dates.format(crtime, "yyyy-mm-dd HH:MM:ss")
    push!(timestamps, (mkt_fips, "end", timestr))
    #print(timestamps, "\n")
    writedlm(outf, timestamps)
    close(outf)
    # Write out temporary files in case the process fails later.
    tout = convert(DataFrame, container);
    names!(tout, colnames)
    tout = tout[ tout[:fipscode].>0, :]
    DataFrames.writetable(pathprograms*"temp_results_$mkt_fips.csv", tout) # no slash.
end


# Add column names to the new data:
#=
output1 = convert(DataFrame, container);

names!(output1, colnames)
output1 = output1[ output1[:fipscode].>0 ,:]

# Append new data to the existing data - there is a stupid problem here that the types don't match exactly.  Make everything Float64:
for el in names(fout1)
  if (typeof(fout1[el]) == DataArrays.DataArray{Int64,1})
#    print("int type ", el, "\n")
    fout1[DataFrames.isna(fout1[el]),:] = -666
  elseif (typeof(fout1[el]) == DataArrays.DataArray{Float64,1})
#    print("float type ", el, "\n")
    fout1[DataFrames.isna(fout1[el]),:] = -666.0
  else
#    print("other type ", typeof(el), "\n")
  end
end

for el in names(fout1)
#  print(el, "\n")
  fout1[el] = convert( Array{Float64,1}, fout1[el])
#  output1[el] = convert( Array{Float64,1}, output1[el])
end

append!(fout1, output1)
fout1 = fout1[fout1[:fipscode].!=-666,:]
writetable(pathprograms*"simulationresults.csv", output1)


=#













###



#=

# Data Types - currently unused.

type Hosp
	fid::Int
	name::AbstractString
	level::Tuple
	fips::Int
	lat::Float64
	long::Float64
	n05::Array
	n515::Array
	n1525::Array
	choices::Array
	probs::WeightVec
end

h1 = Hosp(dataf[1,:fid], dataf[1,:facility], (dataf[1,:act_int], dataf[1,:act_solo]), dataf[1,:fipscode], dataf[1, :v15], dataf[1, :v16],     [dataf[1,:lev105], dataf[1,:lev205], dataf[1,:lev305]], [dataf[1,:lev1515], dataf[1,:lev2515], dataf[1,:lev3515]], [dataf[1,:lev11525], dataf[1,:lev21525], dataf[1,:lev31525]], [1, 2, 3, 4], WeightVec([0.1, 0.1, 0.1, 0.1])     )

type Market
	fips::Int
	lev1::Int
	lev2::Int
	lev3::Int
	config::Array{Hosp}
end


=#
