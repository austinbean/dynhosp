



using DataFrames
using Distributions
using StatsBase
using Plots
plotlyjs()   # call PlotljJS backend to Plots.

#dat = readtable("/Users/austinbean/Google Drive/Simulation Results/dynhospsimresults.csv");
#dat = readtable("/Users/austinbean/Desktop/dynhospsimulationresults2.csv"); # This one is just for testing purposes.  Not a real set of results.

dat = readtable("/Users/austinbean/Desktop/dynhospsimulationresults 10 03 2016 717pm.csv");

dat = convert(Array{Float64,2}, dat)

# Equilibrium and non-equilibrium Medicaid patient revenue:
eq_const = sum( dat[:, 26:31], 2);
neq_const = sum( dat[:,66:72], 2);

# Equilibrium and non-equilibrium WTP, DRG costs, per-patient values:
interimeq_opt = hcat(dat[:, 2:25], dat[:, 33:41]);
interimneq_opt = hcat(dat[:,42:65], dat[:,73:81]);


function objfun(x::Vector;
                scale_fact = 1,
                inp1::Array{Float64,2}=scale_fact*interimeq_opt,
                inp2::Array{Float64,2}=scale_fact*interimneq_opt,
                cons1::Array{Float64,2}=scale_fact*eq_const,
                cons2::Array{Float64,2}=scale_fact*neq_const,
                diffmat::Array{Float64,2}=inp1-inp2)
     """
   this is the BBL objective function
     """
  return sum(min(diffmat*x+eq_const - neq_const, 0).^2)
end

function ObjBehavior(index::Int64, lowerbound::Float64, upperbound::Float64; cnts::Int64 = 1000, argsize::Int64 = 33)
  # takes the objective function above and evaluates it variable-at-a-time between lowerbound and upperbound at cnts points.
  lins = linspace(lowerbound, upperbound, cnts)
  inp = zeros(argsize, cnts)
  outp = zeros(cnts)
  for i = 1:size(lins, 1)
    inp[index, i] = lins[i]
    outp[i] = objfun(inp[:,i])
  end
  return outp, lins
end

function PlotAll(lb::Float64, ub::Float64; argsize::Int64 = 33, cnts::Int64 = 1000)
  outp = zeros(argsize, cnts)
  for i = 1:argsize
    out1, inp1 = ObjBehavior(i, lb, ub)
    outp[i,:] = out1
  #  plot(outp[i,:], linspace(lb, ub, cnts))
  end
  return outp, linspace(lb, ub, cnts)
end

A1, A2 = PlotAll(-100.0, 100.0);
p1 = Plot(plot(A2, A1[1,:]))
p2 = Plot(plot(A2, A1[2,:]))
p3 = plot(A2, A1[3, :])
p4 = plot(A2, A1[4, :])
p5 = plot(A2, A1[5, :])
p6 = plot(A2, A1[6, :])
p7 = plot(A2, A1[7, :])
p8 = plot(A2, A1[8, :])
p9 = plot(A2, A1[9, :])
p10 = plot(A2, A1[10, :])
p11 = plot(A2, A1[11, :])
p12 = plot(A2, A1[12, :])
p13 = plot(A2, A1[13,:])
p14 = plot(A2, A1[14,:])
p15 = plot(A2, A1[15,:])
p16 = plot(A2, A1[16,:])
p17 = plot(A2, A1[17,:])
p18 = plot(A2, A1[18,:])
p19 = plot(A2, A1[19,:])
p20 = plot(A2, A1[20,:])
p21 = plot(A2, A1[21,:])
p22 = plot(A2, A1[22,:])
p23 = plot(A2, A1[23, :])
p24 = plot(A2, A1[24, :])
p25 = plot(A2, A1[25, :])
p26 = plot(A2, A1[26, :])
p27 = plot(A2, A1[27, :])
p28 = plot(A2, A1[28, :])
p29 = plot(A2, A1[29, :])
p30 = plot(A2, A1[30, :])
p31 = plot(A2, A1[31, :])
p32 = plot(A2, A1[32, :])
p33 = plot(A2, A1[33,:])
