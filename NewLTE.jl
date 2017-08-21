# Redoing the LTE.



#TODO - remove DataFrames, sub for readcsv.  Ugh.  
#using DataFrames
using Distributions
using StatsBase
using Plots
plotlyjs()   # call PlotlyJS backend to Plots.

#dat = readtable("/Users/austinbean/Google Drive/Simulation Results/dynhospsimresults.csv");
#dat = readtable("/Users/austinbean/Desktop/dynhospsimulationresults2.csv"); # This one is just for testing purposes.  Not a real set of results.

dat = readcsv("/Users/austinbean/Desktop/dynhosp/Results/longobjective2017-07-21-19-03-37.csv", header = false)

#dat = convert(Array{Float64,2}, dat)

# Equilibrium and non-equilibrium Medicaid patient revenue:
eq_const = sum( dat[:, 26:31], 2);
neq_const = sum( dat[:,66:72], 2);

# Equilibrium and non-equilibrium WTP, DRG costs, per-patient values:
interimeq_opt = hcat(dat[:, 2:25], dat[:, 33:41]);
interimneq_opt = hcat(dat[:,42:65], dat[:,73:81]);

"""
`objfun(x::Vector;scale_fact = 1,inp1::Array{Float64,2}=scale_fact*interimeq_opt,inp2::Array{Float64,2}=scale_fact*interimneq_opt,cons1::Array{Float64,2}=scale_fact*eq_const,cons2::Array{Float64,2}=scale_fact*neq_const,diffmat::Array{Float64,2}=inp1-inp2)`

The BBL objective function.

## Testing ## 

bblobjfun(ones(33), interimeq_opt, interimneq_opt, eq_const, neq_const)

TODO - I think I want this to be MINIMIZED, so should return -1 of this.  
"""
function bblobjfun(x::Vector, inp1::Array{Float64,2}, inp2::Array{Float64,2}, cons1::Array{Float64,2}, cons2::Array{Float64,2})
  # this is the BBL objective function
  nc::Float64 = 1/size(inp1,1)
  return -nc*sum(min.((inp1-inp2)*x+cons1- cons2, 0).^2)
  # this is not the smartest way to do this anyway: why not loop over elements?
end


"""
bbl2(ones(33), interimeq_opt, interimneq_opt, eq_const, neq_const)

Not done yet.  Gives different answer.  But allocation is lower and speed is better.  


I think I want this to be MINIMIZED, so should return -1 of this. 

Ok - time to debug this. 

# test this.
sm = 0.0  
for j = 1:33
  sm += (interimeq_opt[1,j]-interimneq_opt[1,j]) 
end 
println(sm)
println(eq_const[1], " ", neq_const[1])
sm += (eq_const[1] - neq_const[1])

println(-(1/286)*(min(sm,0 )^2) )

-3.3546970628060017e18

"""
function bbl2(x::Vector, inp1::Array{Float64,2}, inp2::Array{Float64,2}, cons1::Array{Float64,2}, cons2::Array{Float64,2})
  nc::Float64 = 1/size(inp1,1)
  sm::Float64 = 0.0
  interim::Float64 = 0.0
  for i = 1:size(inp1,1)                                   # rows of inp1 - single inequalities.  
    interim = 0.0
    for j = 1:size(x,1)                                    # rows of x - parameter values.
      interim += (inp1[i,j].-inp2[i,j]).*x[j]
    end 
    if i == 1 println(interim) end # correct through here.
    if i == 1 println(cons1[i], "  ", cons2[i]) end
    sm += 1/nc*((min(interim+cons1[i]-cons2[i],0))^2)      # params*(eq_opt - neq_opt) + eq_const - neq_const
    if i == 1
      println("sm ", -sm)
      println("at 1: ", -1/nc*(min(interim+cons1[1]-cons2[1],0)^2))
    end 
  end 
  return -sm 
end 




"""
`TestMH`

guess = [100.0, 100.0, 100.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412, 2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ];

testfun(guess) # should be 0.0

"""
function testfun(x::Array{Float64})
  #              [100.0,          100.0,            100.0,          12038.0,           12038,           12038,           66143,           66143,           66143,            19799,            19799,            19799,            4044,            4044,            4044,            6242,            6242,            6242,            1329,            1329,            1329,            412,            412,            412,            2000000,            5000000,            0,           -100000,            2000000,            0,           -200000,           -100000,            0 ],
  return -(((x[1]-100.0)^2)+((x[2]-100.0)^2)+((x[3]-100.0)^2)+((x[4]-12038.0)^2)+((x[5]-12038)^2)+((x[6]-12038)^2)+((x[7]-66143)^2)+((x[8]-66143)^2)+((x[9]-66143)^2)+((x[10]-19799)^2)+((x[11]-19799)^2)+((x[12]-19799)^2)+((x[13]-4044)^2)+((x[14]-4044)^2)+((x[15]-4044)^2)+((x[16]-6242)^2)+((x[17]-6242)^2)+((x[18]-6242)^2)+((x[19]-1329)^2)+((x[20]-1329)^2)+((x[21]-1329)^2)+((x[22]-412)^2)+((x[23]-412)^2)+((x[24]-412)^2)+((x[25]-2000000)^2)+((x[26]-5000000)^2)+((x[27]-0)^2)+((x[28]+100000)^2)+((x[29]-2000000)^2)+((x[30]-0)^2)+((x[31]+200000)^2)+((x[32]+100000)^2)+((x[33]-0)^2))
end 





""" 
`GetModes(x::Array{Float64,2})`

Compute the mode of the values:
"""
function GetModes(x::Array{Float64,2})
  # Round these values and get the modes of the round values
  outp = zeros(size(x,2))
  for i = 1:size(x,2)
    outp[i] = mode(round.(x[:,i], 4)) # round to four digits and take the mode.
  end
  return outp
end



"""
`ResultsPrint(x::Array{Float64,1}, start::Array{Float64,1})`
Take the initial guess and the final value and print them next to the parameter name.
"""
function ResultsPrint(x::Array{Float64,1}, start::Array{Float64,1})
  syms = [:α₁, :α₂, :α₃, :γ¹₅, :γ₅², :γ₅³, :γ₆¹, :γ₆², :γ₆³, :γ₇¹, :γ₇², :γ₇³, :γ₈¹, :γ₈², :γ₈³, :γ₉¹, :γ₉², :γ₉³, :γ₀¹, :γ₀², :γ₀³, :γ₁¹, :γ₁², :γ₁³, :ϕ12, :ϕ13, :ϕ1EX, :ϕ21, :ϕ23, :ϕ2EX, :ϕ31, :ϕ32, :ϕ3EX]
  for el in 1:size(x,1)
    println( syms[el], "  Initial: ", round(start[el],2), "  Final:  ", round(x[el], 2))
  end
end


const nsims = 100 #_000

#drgamt::Array{Float64,1} = [12038.83, 66143.19, 19799.52, 4044.67, 6242.39, 1329.98, 412.04]
guess = [10.0, 1.0, 1.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412, 2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ];
# sim_vals, overcounter, undercounter, accept, tr, param_accept, allvals = MetropolisHastings(guess, nsims; debug = true) # for debugging


sim_vals, overcounter, undercounter, accepted = MetropolisHastings(guess, nsims, testfun; debug = false) # no debugging output.


sim_vals, overcounter, undercounter, accepted = MetropolisHastings(guess, nsims, testfun ; debug = false) # no debugging output.



# Results to return:
println("Fraction Accepted ", accepted/nsims)
println("Count of Numerical Overflow ", overcounter)
println("Count of Underflow ", undercounter)


# NB!! Re-enter "Guess"

ans1 = GetModes(sim_vals);
guess = [10.0, 1.0, 1.0, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412, 2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ];
ResultsPrint(ans1, guess)


# Plotting Results:

histogram(sim_vals[:,1], nbins=20)
histogram(sim_vals[:,2], nbins=20)


# Sim Annealing:
guess = [1.0, 1, 1, 12038.0, 12038, 12038, 66143, 66143, 66143, 19799, 19799, 19799, 4044, 4044, 4044, 6242, 6242, 6242, 1329, 1329, 1329, 412, 412, 412, 2000000, 5000000, 0, -100000, 2000000, 0, -200000, -100000, 0 ];

res1 =  optimize(objfun, guess, method = SimulatedAnnealing(), iterations = 10_000_0, show_trace = true, show_every = 500_000)



#=
# These columns just record what is where.
arr[1] = (alph1 = (beta^T)*outprob*dot(wtp_out[1:7], private[1:7])  )           # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 1
arr[2] = (alph2 = (beta^T)*outprob*dot(wtp_out[8:14], private[8:14])   )        # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 2
arr[3] = (alph3 = (beta^T)*outprob*dot(wtp_out[15:21], private[15:21]) )        # This is WTP over all DRGS * patient vols at corresponding DRG, over all periods at level 3
arr[4] = (gamma_1_385 = (beta^T)*outprob*(private[1]+medicaid[1])   )           # The next lines are patients summed over types.  Costs are treated as the same over Medicaid and privately insured.
arr[5] = (gamma_2_385 = (beta^T)*outprob*(private[8]+medicaid[8]) )             # I don't use any of the named terms - they just keep track.
arr[6] = (gamma_3_385 = (beta^T)*outprob*(private[15]+medicaid[15]) )
arr[7] = (gamma_1_386 = (beta^T)*outprob*(private[2]+medicaid[2]) )
arr[8] = (gamma_2_386 = (beta^T)*outprob*(private[9]+medicaid[9]) )
arr[9] = (gamma_3_386 = (beta^T)*outprob*(private[16]+medicaid[16]) )
arr[10] = (gamma_1_387 = (beta^T)*outprob*(private[3]+medicaid[3]) )
arr[11] = (gamma_2_387 = (beta^T)*outprob*(private[10]+medicaid[10]) )
arr[12] = (gamma_3_387 = (beta^T)*outprob*(private[17]+medicaid[17]) )
arr[13] = (gamma_1_388 = (beta^T)*outprob*(private[4]+medicaid[4]) )
arr[14] = (gamma_2_388 = (beta^T)*outprob*(private[11]+medicaid[11]) )
arr[15] = (gamma_3_388 = (beta^T)*outprob*(private[18]+medicaid[18]) )
arr[16] = (gamma_1_389 = (beta^T)*outprob*(private[5]+medicaid[5]) )
arr[17] = (gamma_2_389 = (beta^T)*outprob*(private[12]+medicaid[12]) )
arr[18] = (gamma_3_389 = (beta^T)*outprob*(private[19]+medicaid[19]) )
arr[19] = (gamma_1_390 = (beta^T)*outprob*(private[6]+medicaid[6]) )
arr[20] = (gamma_2_390 = (beta^T)*outprob*(private[13]+medicaid[13]) )
arr[21] = (gamma_3_390 = (beta^T)*outprob*(private[20]+medicaid[20]) )
arr[22] = (gamma_1_391 = (beta^T)*outprob*(private[7]+medicaid[7]) )
arr[23] = (gamma_2_391 = (beta^T)*outprob*(private[14]+medicaid[14]) )
arr[24] = (gamma_3_391 = (beta^T)*outprob*(private[21]+medicaid[21]) )
arr[25] = (patients385 = (beta^T)*outprob*(medicaid[1]+medicaid[8]+medicaid[15]))    # Count of patients at DRG 385
arr[26] = (patients386 = (beta^T)*outprob*(medicaid[2]+medicaid[9]+medicaid[16]))
arr[27] = (patients387 = (beta^T)*outprob*(medicaid[3]+medicaid[10]+medicaid[17]))
arr[28] = (patients388 = (beta^T)*outprob*(medicaid[4]+medicaid[11]+medicaid[18]))
arr[29] = (patients389 = (beta^T)*outprob*(medicaid[5]+medicaid[12]+medicaid[19]))
arr[30] = (patients390 = (beta^T)*outprob*(medicaid[6]+medicaid[13]+medicaid[20]))
arr[31] = (patients391 = (beta^T)*outprob*(medicaid[7]+medicaid[14]+medicaid[21]))
arr[32:end] = (alltrans = (beta^T)*outprob*transitions)

=#
