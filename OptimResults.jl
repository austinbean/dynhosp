# Take the outcome of main and give it to the optimizer:

using DataFrames
using DataArrays
using Distributions
using Optim
using PyPlot


fout1 = readtable("/Users/austinbean/Desktop/dynhosp/simulationresults.csv")



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

names!(fout1, colnames)

paramsymbs = Array{UTF8String}(0)
for k = 1:3
  for i = 0:25
    push!(paramsymbs, "Θ"*"$k"*"C$i")
  end
end
for y in [1 2 3]
  for z in [1 2 3 "EX"]
    if y != z
      push!(paramsymbs, "ψ"*"$y"*"$z")
    end
  end
end

push!(paramsymbs, "Γ1")
push!(paramsymbs, "Γ2")
push!(paramsymbs, "Γ3")

for i in 1:size(paramsymbs)[1]
  print(colnames[3+i],"  ", paramsymbs[i], "\n")
end





function dfvec(datafr::DataFrame)
  vecvals = Vector{Float64}(0)
  for el in 1:size(datafr)[2]
    push!(vecvals, datafr[el].data[1])
  end
  return vecvals
end

# Delete columns of zeros - this is just for the testing part.  Eventually hopefully all will be filled in.
# This only deletes if both the column for the equilibrium AND non-equilibrium are zero.
# this doesn't work because the size is changing dynamically.
deletdsym = Array{Int64}(0)
deletdcols = Array{Int64}(0)
for col in colnames[4:end-90]
  el = fout1.colindex.lookup[col] #gets column number from dictionary
  nel = el + 90
  if (sum(fout1[colnames[el]]) == 0) & (sum(fout1[colnames[el+90]]) == 0)
    print("Empty Column: ", colnames[el], " ", colnames[el+90], " Symbol:", paramsymbs[el-3] ,"\n")
    push!(deletdcols, el)
    push!(deletdcols, nel)
    push!(deletdsym, el-3)
  end
end
deleteat!(paramsymbs, deletdsym)
delete!(fout1, deletdcols)

# Number of simulations and num
hsims = size(fout1)[1] # number of simulations
ncols = size(fout1)[2] # number of columns with nonzeros (pairs!)

# How many parameters am I trying to estimate?  # of Non-zero column pairs
params = convert(Int, (ncols-3)/2) # don't think this conversion is strictly necessary

# Check resulting:
for el in 1:params
  print(fout1.colindex.names[el+3], "  ", paramsymbs[el], "\n")
end



for x in 1:hsims
  name = parse("function val$x(x::Vector; inp1 = dfvec(fout1[$x, 4:4+params-1]), inp2 = dfvec(fout1[$x, 4+params:end])) return (minimum([sum(x.*(inp1 - inp2)), 0.0]))^2 end")
  eval(name)
end

#=
To evaluate the above:
for x in 1:20
 phrs = parse("val$x(ones(params))")
 print(eval(phrs), "\n")
end
=#
str = ""
for x in 1:hsims
 str = str*"val$x(x) + "
end

function sumval(x::Vector; hsims = 75)
  return (1/hsims)*(val1(x) + val2(x) + val3(x) + val4(x) + val5(x) + val6(x) + val7(x) + val8(x) + val9(x) + val10(x) + val11(x) + val12(x) + val13(x) + val14(x) + val15(x) + val16(x) + val17(x) + val18(x) + val19(x) + val20(x) + val21(x) + val22(x) + val23(x) + val24(x) + val25(x) + val26(x) + val27(x) + val28(x) + val29(x) + val30(x) + val31(x) + val32(x) + val33(x) + val34(x) + val35(x) + val36(x) + val37(x) + val38(x) + val39(x) + val40(x) + val41(x) + val42(x) + val43(x) + val44(x) + val45(x) + val46(x) + val47(x) + val48(x) + val49(x) + val50(x) + val51(x) + val52(x) + val53(x) + val54(x) + val55(x) + val56(x) + val57(x) + val58(x) + val59(x) + val60(x) + val61(x) + val62(x) + val63(x) + val64(x) + val65(x) + val66(x) + val67(x) + val68(x) + val69(x) + val70(x) + val71(x) + val72(x) + val73(x) + val74(x) + val75(x))
end


result = optimize(sumval, ones(params), method = SimulatedAnnealing(), iterations = 5000)

result = optimize(sumval, ones(params), method = SimulatedAnnealing(), iterations = 50000, store_trace = true)

result2 = optimize(sumval, ones(params), method = SimulatedAnnealing(), iterations = 5000, extended_trace = true, store_trace = true);

result = optimize(sumval, 500*ones(params), method = SimulatedAnnealing(), iterations = 50000, store_trace = true, ftol = 1.0e-8)


# Now this will print the parameter name:
for el in 1:convert(Int, (size(fout1.colindex.names)[1]-3)/2)
  print(fout1.colindex.names[el+3], "  ", Optim.minimizer(result)[el], " param: ", paramsymbs[el], "\n")
end

# These are not very informative yet.  
plot(Optim.minimizer(result)[1:4], paramsymbs[1:4])
plot(Optim.minimizer(result)[5:8], paramsymbs[5:8])


##


#=
# Test from the Optim Page::
function f(x::Vector)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
optimize(f, [0.0, 0.0])

Outcomes:

* Algorithm: Nelder-Mead
 * Starting Point: [0.0,0.0]
 * Minimizer: [1.000005438687492,1.0000079372595394]
 * Minimum: 0.000000
 * Iterations: 60
 * Convergence: true
   * |x - x'| < NaN: false
   * |f(x) - f(x')| / |f(x)| < 1.0e-08: true
   * |g(x)| < NaN: false
   * Reached Maximum Number of Iterations: false
 * Objective Function Calls: 115
 * Gradient Calls: 0


=#
