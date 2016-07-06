# Take the outcome of main and give it to the optimizer:

using DataFrames
using Distributions
using Optim
using Gadfly

# fout1 = readtable("/Users/austinbean/Desktop/dynhosp/simulationresults.csv")
#fout1 = readtable("/Users/austinbean/Desktop/temp_results_48209.csv")

fout1 = readtable("/Users/austinbean/Google Drive/Simulation Results/simulationresults.csv");


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

#names!(fout1, colnames)

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
  for el in 1:size(datafr,2)
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
hsims = 500 #size(fout1)[1] # number of simulations
ncols = size(fout1,2) # number of columns with nonzeros (pairs!)

# How many parameters am I trying to estimate?  # of Non-zero column pairs
params = convert(Int, (ncols-3)/2) # don't think this conversion is strictly necessary

# Check resulting:
for el in 1:params
  print(fout1.colindex.names[el+3], "  ", paramsymbs[el], "\n")
end



for x in 1:size(fout1,1) #hsims #how many facilities are we doing?  Not hsims
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
for x in 1:size(fout1,1)
 str = str*"val$x(x) + "
end

function sumval(x::Vector; hsims = 500)
  return (1/hsims)*(val1(x) + val2(x) + val3(x) + val4(x) + val5(x) + val6(x) + val7(x) + val8(x) + val9(x) + val10(x) + val11(x) + val12(x) + val13(x) + val14(x) + val15(x) + val16(x) + val17(x) + val18(x) + val19(x) + val20(x))
end


result = optimize(sumval, ones(params), method = SimulatedAnnealing(), iterations = 5000)

result = optimize(sumval, ones(params), method = SimulatedAnnealing(), iterations = 50000, store_trace = true)

result2 = optimize(sumval, ones(params), method = SimulatedAnnealing(), iterations = 5000, extended_trace = true, store_trace = true);

result = optimize(sumval, 500*ones(params), method = SimulatedAnnealing(), iterations = 50000, store_trace = true)

result3 = optimize(sumval, 700*ones(params), method = SimulatedAnnealing(), iterations = 50000, store_trace = true)


# Now this will print the parameter name:
for el in 1:convert(Int, (size(fout1.colindex.names)[1]-3)/2)
  print(fout1.colindex.names[el+3], "  ", Optim.minimizer(result)[el], " param: ", paramsymbs[el], " symbol number: ", el, "\n")
end

# These are not very informative yet.

# Test labeling:
x1 = paramsymbs[1:20]
p1 = plot(x=x1, y=Optim.minimizer(result3)[1:20], Guide.xticks(ticks=collect(1:20)), Guide.xlabel("Profits Per Patient at Level 1 \n Given Competitors C#"), Guide.ylabel("Dollars"), Geom.point, Geom.line ) # the collect[1:6] is the NUMBER of parameters, not their indices
draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/lev1rev.png", 12cm, 6cm), p1)

x2 = paramsymbs[7:10]
p2 = plot(x=x2, y=Optim.minimizer(result)[7:10], Guide.xticks(ticks=collect(1:4)), Guide.xlabel("Profits Per Patient at Level 2 \n Given Competitors C#"), Guide.ylabel("Dollars"), Geom.point, Geom.line  )
draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/lev2rev.png", 12cm, 6cm), p2)

x3 = paramsymbs[11:13]
p3 = plot(x=x3, y=Optim.minimizer(result)[11:13], Guide.xticks(ticks=collect(1:3)), Guide.xlabel("Profits Per Patient at Level 3 \n Given Competitors C#"), Guide.ylabel("Dollars"), Geom.point, Geom.line  )
draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/lev3rev.png", 12cm, 6cm), p3)

x4 = collect(1:6)
p4 = plot(
layer(x=x4, y=Optim.minimizer(result)[1:6], Geom.point, Geom.line, Theme(default_color=color("green")) ),
layer(x=collect(1:4), y=Optim.minimizer(result)[7:10], Geom.point, Geom.line, Theme(default_color=color("purple"))),
layer(x=collect(1:3), y=Optim.minimizer(result)[11:13], Geom.point, Geom.line, Theme(default_color=color("red"))),
Guide.xlabel("Number of Competitors"),
Guide.ylabel("Dollars"),
Guide.title("Profit per Patient at Three Levels as a Function of Number of Competitors"),
Guide.manual_color_key("Levels", ["Level 1", "Level 2", "Level 3"], ["green", "purple", "red"]))
draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/alllevs.png", 12cm, 6cm), p4)


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
