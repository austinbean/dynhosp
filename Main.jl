# Main

#=
Julia Version 0.4.5
Commit 2ac304d (2016-03-18 00:58 UTC)
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




### Collect Basic Information ###

# locates starting and ending points of all separate fipscode indices



# Next one stores all the years appearing in each fips code
#yearins = [ [x; findfirst(data[:,fipscodeloc], x); findlast(data[:,fipscodeloc], x ); unique( data[findfirst(data[:,fipscodeloc], x):findlast(data[:,fipscodeloc], x ) , yearloc]  ) ] for x in unique(data[:,fipscodeloc])  ];

yearins = [ x   for x in unique(data[ data[:,yearloc].==2005 , fipscodeloc])]

monopoly = Array{Int64}(0)
duopoly = Array{Int64}(0)
triopoly = Array{Int64}(0)
tetrapoly = Array{Int64}(0)
nopoly = Array{Int64}(0)
all = Array{Int64}(0)

for el in yearins
#  unqfids = [x for x in unique(data[el[2]:el[3],fidloc])]
  unqfids = [x for x in unique(data[ (data[:,fipscodeloc].==el[1])&(data[:,yearloc].==2005) ,fidloc])]
  if size(unqfids,1) == 1
  #  print("Fipscode Monopoly: ", unique(dataf[el[2]:el[3], :fipscode]), "\n")
    println("Fipscode Monopoly: ", el[1])
    push!(monopoly, el[1])
  elseif size(unqfids,1) == 2
    println("Fipscode Duopoly: ", el[1])
    push!(duopoly, el[1])
  elseif size(unqfids,1) == 3
    println("Fipscode Triopoly: ", el[1])
    push!(triopoly, el[1])
  elseif size(unqfids,1) == 4
    println("Fipscode Tetrapoly: ", el[1])
    push!(tetrapoly, el[1])
  elseif size(unqfids,1) > 4
    println("Fipscode N-opoly: ", el[1])
    println("Fipscode Hospitals: ", size(unqfids, 1))
    push!(nopoly, el[1])
  end
  push!(all, size(unqfids,1))
end


container = zeros(1, ((2*78+12)*2)+3) #must agree with outp in ParMainFun

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
  for name in ["PI", "MED"]
    for j = 1:3
      for k = 0:25
        push!(colnames, parse("$name""$elem"*"Lev$j"*"Comp$k"))
      end
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

for y in 1:size(duopoly,1)
    mkt_fips = duopoly[y][1]
    crtime = now()
    timestr = Dates.format(crtime, "yyyy-mm-dd HH:MM:ss")
    push!(timestamps, (mkt_fips, "begin", timestr))
    if !(in(mkt_fips, donefips)) & (mkt_fips != 48201) & (mkt_fips != 0) # this is going to do new fipscodes only
      print("Market FIPS Code ", mkt_fips, "\n")
      	for year in [ 2005 ]   #yearins[y][4:end] # can do all years or several.
          fids =  sort!(convert(Array{Int64}, unique(data[(data[:,fipscodeloc].==mkt_fips)&(data[:, yearloc].==year),fidloc])))
          # This will parallelize the computation across Monte Carlo sims.
          mcres = @parallel (+) for i = 1:500
    #                println(i)
                    ParMainfun(data, pinsured, privatedemandmodelparameters, pmedicaid, medicaiddemandmodelparameters, mkt_fips, year,  fids; npers = 50)
                  end
          # the addition map adds together years, fipscodes and fids.  Divide by nsims to recover.
          mcres[:, 1:3] = mcres[:,1:3]/500
          container = [container; mcres]
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
output1 = convert(DataFrame, container);
names!(output1, colnames)
output1 = output1[ output1[:fipscode].>0 ,:]

# Append new data to the existing data - there is a stupid problem here that the types don't match exactly.  Make everything Float64:
for el in names(fout1)
  if (typeof(fout1[el]) == DataArrays.DataArray{Int64,1})
    fout1[DataFrames.isna(fout1[el]),:] = -666
  elseif (typeof(fout1[el]) == DataArrays.DataArray{Float64,1})
    fout1[DataFrames.isna(fout1[el]),:] = -666.0
  else
#    print("other type ", typeof(el), "\n")
  end
end

for el in names(fout1)
  fout1[el] = convert( Array{Float64,1}, fout1[el])
end

append!(fout1, output1)
fout1 = fout1[fout1[:fipscode].!=-666,:]
writetable(pathprograms*"simulationresults.csv", output1)
















###



#=

allindices = [ (x, findfirst(data[:,fipscodeloc], x), findlast(data[:,fipscodeloc], x)) for x in unique(data[:,fipscodeloc]) ];


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
=#
