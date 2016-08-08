# Take the outcome of main and give it to the optimizer:
# include("/Users/austinbean/Desktop/dynhosp/OptimResults.jl")


using DataFrames
using Distributions
using Optim
using Gadfly
using ForwardDiff

# fout1 = readtable("/Users/austinbean/Desktop/dynhosp/simulationresults.csv")
#fout1 = readtable("/Users/austinbean/Desktop/temp_results_48209.csv")
function mainopt()
  input1 = readtable("/Users/austinbean/Google Drive/Simulation Results/combinedresults.csv");

  # This keeps the privately insured patient counts:
  # The division between first and second half is around 172.
    interim_private_eq = convert(Matrix, hcat(input1[:,1:81], input1[:,160:171]));
    interim_private_neq = convert(Matrix, hcat(input1[:,1:3], input1[:,172:249], input1[:,328:339]));
  # This keeps the Medicaid patient counts:
    interim_medicaid_eq = hcat(input1[:,1:3], input1[:,82:159]);
    interim_medicaid_neq = hcat(input1[:,1:3], input1[:,250:327]);
  # This keeps the equilibrium sim medicaid patient counts by level - interim_medicaid summed across competitor numbers
    eq_medicaid_lev1 = hcat(convert(Matrix, interim_medicaid_eq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_eq[:,4:29]),2));
    eq_medicaid_lev2 = hcat(convert(Matrix, interim_medicaid_eq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_eq[:,30:55]),2));
    eq_medicaid_lev3 = hcat(convert(Matrix, interim_medicaid_eq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_eq[:,56:81]),2));
  # This keeps the non-equilibrium sim medicaid counts by level
    neq_medicaid_lev1 = hcat(convert(Matrix, interim_medicaid_neq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_neq[:,4:29]),2));
    neq_medicaid_lev2 = hcat(convert(Matrix, interim_medicaid_neq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_neq[:,30:55]),2));
    neq_medicaid_lev3 = hcat(convert(Matrix, interim_medicaid_neq[:,1:3]) ,sum(convert(Matrix, interim_medicaid_neq[:,56:81]),2));


  # The form of the output is the same as the input, except at each level the 26 medicaid columns have been summed into one
  # This creates TWO matrices, one with all of the equilibrium results, the other with all of the non-equilibrium results.
  fout11 = hcat( interim_private_eq[:,1:29], eq_medicaid_lev1[:,4], interim_private_eq[:,30:55], eq_medicaid_lev2[:,4], interim_private_eq[:,56:81], eq_medicaid_lev3[:,4], interim_private_eq[:,82:end])
  fout12 = hcat(interim_private_neq[:,1:29], neq_medicaid_lev1[:,4], interim_private_neq[:,30:55], neq_medicaid_lev2[:,4], interim_private_neq[:,56:81], neq_medicaid_lev3[:,4], interim_private_neq[:,82:end])

  varcolnames = Array{Symbol}(:0)
  push!(varcolnames, :fipscode)
  push!(varcolnames, :fid)
  push!(varcolnames, :year)
  for elem in ["EQ", "NEQ"]
    for name in ["PI"]
      for j = 1:3
        for k = 0:25
          push!(varcolnames, parse("$name"*"$elem"*"Lev$j"*"Comp$k"))
        end
          push!(varcolnames, parse("$elem"*"Medicaid"*"Lev$j"))
      end
    end
    for x in [1 2 3]
      for y in [1 2 3 "EX"]
        if x != y
          push!(varcolnames, parse("$elem"*"Trans$x$y"))
        end
      end
    end
    push!(varcolnames, parse("$elem"*"Enter1"))
    push!(varcolnames, parse("$elem"*"Enter2"))
    push!(varcolnames, parse("$elem"*"Enter3"))
  end

  #names!(fout1, colnames)

  paramsymbs = Array{UTF8String}(0)
  for k = 1:3
    for i = 0:25
      push!(paramsymbs, "Θ"*"$k"*"C$i")
    end
    push!(paramsymbs, "μ"*"$k")
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


  #=

  # This just checks that varcolnames and parasymbs are correct

  for i in 1:size(paramsymbs,1)
    println(varcolnames[3+i],"  ", paramsymbs[i])
  end

  # Delete columns of zeros - this is just for the testing part.  Eventually hopefully all will be filled in.
  # This only deletes if both the column for the equilibrium AND non-equilibrium are zero.
  # this doesn't work because the size is changing dynamically.
  # This will check quickly if there are pairs of columns which are all zeros.

  =#

  deletdsym = Array{Int64}(0) # vector of ints, for symbol indices (no "push!" method for arrays of symbols?)
  deletdcols = Array{Int64,1}(0) # vector of ints, for column indices - can index columns with this.  Start with ALL ints and remove these.
  pcop = deepcopy(paramsymbs)
  for col in 4:size(fout11,2)
    el = varcolnames[col] # these two lines give column names
    nel = varcolnames[col+93]
    if (sum(fout11[:,col]) == 0) & (sum(fout12[:,col]) == 0)
  #    println("Empty Column: ", col, "  ",sum(fout11[:,col]), "  ", el, "   ",sum(fout12[:,col]), "   ", nel )
    #  fout1 = fout1[:, (1:size(fout1,2).!=col)&(1:size(fout1,2).!=col+93)]
      push!(deletdcols, col)
      push!(deletdsym, col)
    else
    #  println(sum(fout11[:,col]), "  ",sum(fout12[:,col]) )
    end
  end

  deleteat!(paramsymbs, [x-3 for x in deletdsym]) # remember that the matrices include a 3-vector of identifiers.
  deleteat!(varcolnames, deletdsym) # removes the names of the colums - useful after the optimization to print the vals
  # Deleting columns:
  nonzers = setdiff(collect(1:size(fout11,2)), deletdcols) # keep only the columns not in the set deletdcols
  fout11 = fout11[:,nonzers]
  fout12 = fout12[:,nonzers]


  #=

  # This prints a list of the remaining coefficients, plus their indices in paramsymbs.

  for el in 1:size(pcop, 1)
    index = findfirst(paramsymbs, pcop[el])
    if index != 0
      println(index, "   ", paramsymbs[index], "   ", pcop[el])
    else
      println("*****   ", pcop[el])
    end
  end


  =#



  # Drop identifiers:

  eq_opt = convert(Array{Float64, 2}, fout11[:,4:end]);
  neq_opt = convert(Array{Float64, 2}, fout12[:, 4:end]);
  opt = eq_opt - neq_opt;

  # The next function returns a vector of length equal to number of sims.
  function objfun(x::Vector; inp1::Array{Float64,2}=eq_opt, inp2::Array{Float64,2}=neq_opt)
     sum((min(inp1*x - inp2*x, 0)).^2)
  end

  # This will compute the gradient - but check to make sure it's doing what you think.
  # g1 = ForwardDiff.gradient(objfun, ones(size(eq_opt,2)), Chunk{10}())
  # # This one also works, actually.
  # h1 = ForwardDiff.hessian(objfun, ones(size(eq_opt,2)), Chunk{10}())

  # But it is also pretty straightforward to just code the gradient and hessian by hand -
  # that would handle any conceivable function evaluation.  Maybe that is smarter.

  # Number of simulations and num
  const hsims = 500 #size(fout1)[1] # number of simulations
  const ncols = size(eq_opt,2) # number of columns with nonzeros (pairs!)
  const params = convert(Int, ncols) # don't think this conversion is strictly necessary


  #
  # result = optimize(objfun_2, ones(params), method = SimulatedAnnealing(), iterations = 50)
  # result = optimize(objfun_2, ones(params), method = SimulatedAnnealing(), iterations = 50000, store_trace = true)
  # result2 = optimize(objfun_2, ones(params), method = SimulatedAnnealing(), iterations = 5000, extended_trace = true, store_trace = true);
  # result = optimize(objfun_2, 500*ones(params), method = SimulatedAnnealing(), iterations = 50000, store_trace = true)
  # This runs very quickly, even with 500,000 evaluations.
  # store_trace
  const starting = 1000
  result3 = optimize(objfun, starting*ones(params), method = SimulatedAnnealing(), iterations = 1_000_000, show_trace = true, show_every = 100000)


  # Now this will print the parameter name:
  for el in 1:size(Optim.minimizer(result3), 1)
    print(varcolnames[el+3], "  ", Optim.minimizer(result3)[el], " param: ", paramsymbs[el], " symbol number: ", el, "\n")
  end

      return starting, paramsymbs, Optim.minimizer(result3)
end #of mainopt() function

starting, paramsymbs, results = mainopt()


  const α₂ = 0.07
  const α₃ = 0.13

  x1 = paramsymbs[1:19]
  x2 = paramsymbs[21:40]
  x3 = paramsymbs[42:61]

  p1 = plot(
            layer(x=x1, y=results[1:19], Geom.point, Geom.line, Theme(default_color=colorant"green")),
            layer(x=x2, y=results[21:40], Geom.point, Geom.line, Theme(default_color=colorant"purple")),
            layer(x=x3, y=results[42:61], Geom.point, Geom.line, Theme(default_color=colorant"red")),
            layer(x=x1, y=starting*ones(size(results[1:19])), Geom.line, Theme(default_color=colorant"black")),
            Guide.xlabel("Number of Competitors"), Guide.ylabel("Dollars"), Guide.manual_color_key("Levels", ["Level 1", "Level 2", "Level 3", "Initial Value"], ["green", "purple", "red", "black"])
            )


x2prime = paramsymbs[21:39]
x3prime = paramsymbs[42:60]
xlab = collect(1:19)
            p2 = plot(
                      layer(x=x1, y=results[1:19], Geom.point, Geom.line, Theme(default_color=colorant"green")),
                      layer(x=xlab, y=results[1:19]+α₂*results[21:39], Geom.point, Geom.line, Theme(default_color=colorant"purple")),
                      layer(x=xlab, y=results[1:19]+α₃*results[42:60], Geom.point, Geom.line, Theme(default_color=colorant"red")),
                      layer(x=xlab, y=starting*ones(size(results[1:19])), Geom.line, Theme(default_color=colorant"black")),
                      Guide.xlabel("Number of Competitors"), Guide.ylabel("Dollars"), Guide.manual_color_key("Levels", ["Level 1", "Level 2", "Level 3", "Initial Value"], ["green", "purple", "red", "black"])
                      )

                      draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/combinedscaledres.png", 18cm, 9cm), p2)


# These are not very informative yet.
#
# # Test labeling:
# x1 = paramsymbs[1:19]
# p1 = plot(x=x1, y=Optim.minimizer(result3)[1:19], Guide.xticks(ticks=collect(1:19)), Guide.xlabel("Profits Per Patient at Level 1 \n Given Competitors C#"), Guide.ylabel("Dollars"), Geom.point, Geom.line ) # the collect[1:6] is the NUMBER of parameters, not their indices
# #draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/lev1rev.png", 12cm, 6cm), p1)
#
# x2 = paramsymbs[21:40]
# p2 = plot(x=x2, y=Optim.minimizer(result3)[21:40], Guide.xticks(ticks=collect(1:19)), Guide.xlabel("Profits Per Patient at Level 2 \n Given Competitors C#"), Guide.ylabel("Dollars"), Geom.point, Geom.line  )
# #draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/lev2rev.png", 12cm, 6cm), p2)
#
# x3 = paramsymbs[42:61]
# p3 = plot(x=x3, y=Optim.minimizer(result3)[42:61], Guide.xticks(ticks=collect(1:19)), Guide.xlabel("Profits Per Patient at Level 3 \n Given Competitors C#"), Guide.ylabel("Dollars"), Geom.point, Geom.line  )
# #draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/lev3rev.png", 12cm, 6cm), p3)
#
# x4 = collect(1:26)
# p4 = plot(
# layer(x=x4, y=Optim.minimizer(result)[1:26], Geom.point, Geom.line, Theme(default_color=colorant"green") ),
# layer(x=x4, y=Optim.minimizer(result)[28:53], Geom.point, Geom.line, Theme(default_color=colorant"purple")),
# layer(x=x4, y=Optim.minimizer(result)[55:80], Geom.point, Geom.line, Theme(default_color=colorant"red")),
# Guide.xlabel("Number of Competitors"),
# Guide.ylabel("Dollars"),
# Guide.title("Profit per Patient at Three Levels as a Function of Number of Competitors"),
# Guide.manual_color_key("Levels", ["Level 1", "Level 2", "Level 3"], ["green", "purple", "red"]))
# #draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/alllevs.png", 12cm, 6cm), p4)
#
# # Medicaid Revenue:
#
# x5 = collect(1:3)
# medvals = [Optim.minimizer(result)[27] Optim.minimizer(result)[54] Optim.minimizer(result)[81] ]
# p5 = plot(x=x5, y=medvals, Geom.point, Geom.line, Guide.xlabel("Facility Level"), Guide.ylabel("Dollars per Patient"), Guide.title("Per Patient Profit - Medicaid"))
#
# # Expected Revenue - Level 1  + αᵢ × Level i
# α₂ = 0.07, α₃ = 0.13,
#
# x4 = collect(1:26)
# p4 = plot(
# layer(x=x4, y=Optim.minimizer(result)[1:26] , Geom.point, Geom.line, Theme(default_color=colorant"green") ),
# layer(x=x4, y=(Optim.minimizer(result)[1:26]+α₂*Optim.minimizer(result)[28:53]), Geom.point, Geom.line, Theme(default_color=colorant"purple")),
# layer(x=x4, y=(y=Optim.minimizer(result)[1:26]+α₃*Optim.minimizer(result)[55:80]), Geom.point, Geom.line, Theme(default_color=colorant"red")),
# Guide.xlabel("Number of Competitors"),
# Guide.ylabel("Dollars"),
# Guide.title("Profit per Patient at Three Levels as a Function of Number of Competitors"),
# Guide.manual_color_key("Levels", ["Level 1", "Level 2", "Level 3"], ["green", "purple", "red"]))
# draw(PNG("/Users/austinbean/Google Drive/Current Projects/!Job Market/!Job Market Paper/alllevs.png", 12cm, 6cm), p4)
#
#
#
#



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


#=

# These functions mapped a dataframe into two vectors.  Since I am using vectors
# directly via the matrix, I don't need to do that.

# This is used below in function definitions
# function dfvec(datafr::DataFrame)
#   vecvals = Vector{Float64}(0)
#   for el in 1:size(datafr,2)
#     push!(vecvals, datafr[el].data[1])
#   end
#   return vecvals
# end
# function dfvec(datafr::Array{Real,2})
#   vecvals = Vector{Float64}(0)
#   for el in 1:size(datafr,2)
#     push!(vecvals, datafr[el].data[1])
#   end
#   return vecvals
# end
#
# # older set of column names:
#
# colnames = Array{Symbol}(:0)
# push!(colnames, :fipscode)
# push!(colnames, :fid)
# push!(colnames, :year)
# for elem in ["EQ", "NEQ"]
#   for name in ["PI", "MED"]
#     for j = 1:3
#       for k = 0:25
#         push!(colnames, parse("$name""$elem"*"Lev$j"*"Comp$k"))
#       end
#     end
#   end
#   for x in [1 2 3]
#     for y in [1 2 3 "EX"]
#       if x != y
#         push!(colnames, parse("$elem"*"Trans$x$y"))
#       end
#     end
#   end
#   push!(colnames, parse("$elem"*"Enter1"))
#   push!(colnames, parse("$elem"*"Enter2"))
#   push!(colnames, parse("$elem"*"Enter3"))
# end
#
# # Relics from an earlier version when I treated each row separately (very, very slow)
#
#
# # How many parameters am I trying to estimate?  # of Non-zero column pairs
# params = convert(Int, (ncols-3)/2) # don't think this conversion is strictly necessary
#
# fout1 = fout1[:,4:end]
# for x in 1:size(fout1,1)
#   name = parse("function val$x(x::Vector; inp1 = fout1[$x, 1:params-1], inp2 = fout1[$x, params:end]) return (minimum([sum(x.*(inp1 - inp2)), 0.0]))^2 end")
#   eval(name)
# end
#
#
# #To evaluate the above:
# for x in 1:20
#  phrs = parse("val$x(ones(params))")
#  print(eval(phrs), "\n")
# end
#
# str = ""
# for x in 1:size(fout1,1)
#  str = str*"val$x(x) + "
# end
#
# function sumval(x::Vector; hsims = 500)
#   return (1/hsims)*(val1(x) + val2(x) + val3(x) + val4(x) + val5(x) + val6(x) + val7(x) + val8(x) + val9(x) + val10(x) + val11(x) + val12(x) + val13(x) + val14(x) + val15(x) + val16(x) + val17(x) + val18(x) + val19(x) + val20(x) + val21(x) + val22(x) + val23(x) + val24(x) + val25(x) + val26(x) + val27(x) + val28(x) + val29(x) + val30(x) + val31(x) + val32(x) + val33(x) + val34(x) + val35(x) + val36(x) + val37(x) + val38(x) + val39(x) + val40(x) + val41(x) + val42(x) + val43(x) + val44(x) + val45(x) + val46(x) + val47(x) + val48(x) + val49(x) + val50(x) + val51(x) + val52(x) + val53(x) + val54(x) + val55(x) + val56(x) + val57(x) + val58(x) + val59(x) + val60(x) + val61(x) + val62(x) + val63(x) + val64(x) + val65(x) + val66(x) + val67(x) + val68(x) + val69(x) + val70(x) + val71(x) + val72(x) + val73(x) + val74(x) + val75(x) + val76(x) + val77(x) + val78(x) + val79(x) + val80(x) + val81(x) + val82(x) + val83(x) + val84(x) + val85(x) + val86(x) + val87(x) + val88(x) + val89(x) + val90(x) + val91(x) + val92(x) + val93(x) + val94(x) + val95(x) + val96(x) + val97(x) + val98(x) + val99(x) + val100(x) + val101(x) + val102(x) + val103(x) + val104(x) + val105(x) + val106(x) + val107(x) + val108(x) + val109(x) + val110(x) + val111(x) + val112(x) + val113(x) + val114(x) + val115(x) + val116(x) + val117(x) + val118(x) + val119(x) + val120(x) + val121(x) + val122(x) + val123(x) + val124(x) + val125(x) + val126(x) + val127(x) + val128(x) + val129(x) + val130(x) + val131(x) + val132(x) + val133(x) + val134(x) + val135(x) + val136(x) + val137(x) + val138(x) + val139(x) + val140(x) + val141(x) + val142(x) + val143(x) + val144(x) + val145(x) + val146(x) + val147(x) + val148(x) + val149(x) + val150(x) + val151(x) + val152(x) + val153(x) + val154(x) + val155(x) + val156(x) + val157(x) + val158(x) + val159(x) + val160(x) + val161(x) + val162(x) + val163(x) + val164(x) + val165(x) + val166(x) + val167(x) + val168(x) + val169(x) + val170(x) + val171(x) + val172(x) + val173(x) + val174(x) + val175(x) + val176(x) + val177(x) + val178(x) + val179(x) + val180(x) + val181(x) + val182(x) + val183(x) + val184(x) + val185(x) + val186(x) + val187(x) + val188(x) + val189(x) + val190(x) + val191(x) + val192(x) + val193(x) + val194(x) + val195(x) + val196(x) + val197(x) + val198(x) + val199(x) + val200(x) + val201(x) + val202(x) + val203(x) + val204(x) + val205(x) + val206(x) + val207(x) + val208(x) + val209(x) + val210(x))
# end

=#
