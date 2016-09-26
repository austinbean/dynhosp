#-----------------------------------------------------------------------------
# load dependencies
#-----------------------------------------------------------------------------

# load packages
using HDF5, JLD           # package to save multidimensional arrays in julia native format
#using PyPlot=#
using Distributions
using Grid
using NLopt
using ArrayViews
using GaussQuadrature
using JuMP
using Autoreload

# set base directory
Base.cd("/home/tobias/Dropbox/Taxi/")

require("programs/julia/mpec/GaussQuadrature.jl")


#Pkg.pin("JuMP",v"0.7.1")



global bootstrap = 0;
# load user-defined functions
arequire("programs/julia/mpec/functions.jl")
arequire("programs/julia/mpec/packages_and_data.jl")




function estparam(active_stop,active_nonstop,inactive_start,inactive_nonstart,breakm,pi)

	m = Model()

	@defVar(m,  0.0001 <= p_est[1:hours,1:shiftl-1,1:num_types,1:nwage] <= 0.9999);                             
	@defVar(m,  0.0001 <= q_est[1:hours,1:num_types] <= 0.9999);
	@defVar(m, -200 <= EV[1:hours,1:shiftl-1,1:num_types,1:nwage] <= 1400);
	@defVar(m, -200 <= EVA[1:hours,1:shiftl-1,1:num_types] <= 1400);
	@defVar(m,0 <= λ_1 <= 20);                                                                    # λ_0 denotes the parts coefficients of the cost function
	@defVar(m,0 <= λ_2 <= 20);                                                                    # λ_1 denotes the parts coefficients of the cost function
	@defVar(m,0 <= λ_0[1:length(unique(pick_o))] <= 250);                                                        # λ_0 denotes the parts coefficients of the cost function
	@defVar(m,-120 <= f[1:2] <= 400);
	@defVar(m, 0.1 <= σ <= 80);                                                                          # variance of the T1EV starting distribution
	@defVar(m, 0.1 <= σ_s <= 80);                                                                         # variance of the T1EV stopping distribution
	@defVar(m,-300 <= μ[1:2] <= 300);                                                                       # the parameters for the value of the outside option

	setValue(λ_1, .01);
	setValue(λ_2, .4);
	setValue(μ[1], 160.0);
	setValue(μ[2], 160.0);


	setValue(σ_s, 22.4);
	setValue(σ, 22.0);


	for h in [unique(pick_o)]
		setValue(λ_0[h], 21.0)	
	end;

	for h in [unique(pick_fine)]
	    setValue(f[h], 60.0)  
	end;

	for h=1:hours,sl=1:shiftl-1, c=1:num_types, a=1:nwage                                             # defining the stopping probabilties
	    setValue(p_est[h,sl,c,a], 0.03*sl)  
	end;

	for h=1:hours,sl=1:shiftl-1, c=1:num_types, a=1:nwage                                           # expected value functions
		setValue(EV[h,sl,c,a], 14)	
	end;

	for h=1:hours,sl=1:shiftl-1, c=1:num_types                                       # expected value functions, integrated out
		setValue(EVA[h,sl,c], 25)	
	end;

	for h=1:hours, c=1:num_types                                                         # starting probabilities
	    setValue(q_est[h,c], 0.6) 
	end;


	# set objective function
	@setNLObjective(m, Max, sum{(1/10000)*(active_stop[h,sl,c,a]*log(p_est[h,sl,c,a])+active_nonstop[h,sl,c,a]*log(1-p_est[h,sl,c,a])),h=1:hours,sl=1:shiftl-1,c=1:num_types, a=1:nwage}+sum{(1/10000)*(inactive_start[h,c]*log(q_est[h,c])+inactive_nonstart[h,c]*log(1-q_est[h,c])),h=1:hours,c=1:num_types})

	# setting constraints:
	# value-constraint, h < 24, sl < smax-1
	for h=1:hours-1,sl=1:shiftl-2, c=1:num_types, a=1:nwage
	    @addNLConstraint(m, EV[h,sl,c,a]==σ*log(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[h+1,sl+1,c])/σ)))
	end

	# value-constraint, h < 24, sl = smax-1
	for h=1:hours,sl=shiftl-1, c=1:num_types, a=1:nwage
	    @addNLConstraint(m, EV[h,sl,c,a]==σ*log(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]])/σ)))
	end

	# value-constraint, h = 24, sl < smax-1
	for h=hours,sl=1:shiftl-2, c=1:num_types, a=1:nwage
	    @addNLConstraint(m, EV[h,sl,c,a]==σ*log(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[1,sl+1,c])/σ)))
	end

	# stopping-prob, h < 24, sl < smax-1
	for h=1:hours-1,sl=1:shiftl-2, c=1:num_types, a=1:nwage
	    @addNLConstraint(m, p_est[h,sl,c,a]==exp(1/σ)/(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[h+1,sl+1,c])/σ)))
	end

	# stopping-prob, h < 24, sl = smax-1
	for h=1:hours,sl=shiftl-1, c=1:num_types, a=1:nwage
	    @addNLConstraint(m, p_est[h,sl,c,a]==exp(1/σ)/(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]])/σ)))
	end

	# stopping-prob, h = 24, sl < smax-1
	for h=hours,sl=1:shiftl-2, c=1:num_types, a=1:nwage
	    @addNLConstraint(m, p_est[h,sl,c,a]==exp(1/σ)/(exp(1/σ)+exp((pi[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-constr_mat[c,h]*f[pick_fine[h]]+EVA[1,sl+1,c])/σ)))
	end

	# starting probabilities h < 24
	for h=1:hours-1, c=1:num_types
	    @addNLConstraint(m, q_est[h,c]==exp((EVA[h+1,1,c]-lease[h+1])/σ_s)/(exp((μ[work_time[h+1]])/σ_s)+exp((EVA[h+1,1,c]-lease[h+1])/σ_s)))
	end;

	# starting probabilities h == 24
	for h=hours, c=1:num_types
	    @addNLConstraint(m, q_est[h,c]==exp((EVA[1,1,c]-lease[1])/σ_s)/(exp((μ[work_time[1]])/σ_s)+exp((EVA[1,1,c]-lease[1])/σ_s)))
	end;

	# integrating out value function with respect to expected wage
	for h=1:hours,sl=1:shiftl-1, c=1:num_types
	    @addNLConstraint(m, EVA[h,sl,c] == sum{pw[a]*EV[h,sl,c,a],a=1:nwage})
	end


	# initiates the model solver
	solve(m)

	mat_p=getValue(p_est)[:,:,:,:];
	mat_EV=getValue(EV)[:,:,:,:];
	mat_q=getValue(q_est)[:,:];
	param=[getValue(λ_0),getValue(λ_1),getValue(λ_2),getValue(σ),getValue(σ_s),getValue(f),getValue(μ)]

	gc()
	return mat_p,mat_q,mat_EV,param 
end;	



#---------------------------------------------------------------------------------------------------------
# for a given set of paremeters this function computes the value function through value function iteration
#---------------------------------------------------------------------------------------------------------


function compute_valuefunc(λ_0,λ_1,λ_2,σ,σ_s,f,μ,pie,breakm::Array,conm = constr_mat)

	# define variables that we need for iteration
	tol = 10e-9
	diff = 3*tol
	num_types = size(conm)[1];
	NEV= zeros(hours,shiftl-1,num_types,nwage);
	p_est = zeros(hours,shiftl-1,num_types,nwage);
	q_est = zeros(hours,num_types);
	

	# set initial conditions for last iteration
	for sl=shiftl-1,h=1:hours, c=1:num_types, a = 1:nwage 
	    if h < hours 
	        nh = h+1
	    else 
	        nh = 1
	    end;  

	    ind = sl == shiftl - 1; 
	    nsl = min(sl + 1,shiftl-1);
	    NEV[h,sl,c,a]=σ*log(exp(1/σ)+exp((pie[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-conm[c,h]*f[pick_fine[h]])/σ))
	end;

	# fill up from the back
	for sl in [shiftl-2:-1:1],h=1:hours, c=1:num_types, a = 1:nwage 

	    if h < hours 
	        nh = h+1
	    else 
	        nh = 1
	    end;  

	    ind = sl == shiftl - 1; 
	    nsl = min(sl + 1,shiftl-1);
	    
	    ANEV = 0.0
	    for ia = 1:nwage
	    	ANEV += pw[ia]*NEV[nh,nsl,c,ia]
	    end;
	    #print(ANEV,"\n")	

	    NEV[h,sl,c,a]=σ*log(exp(1/σ)+exp((pie[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-conm[c,h]*f[pick_fine[h]]+ANEV)/σ))
	end;

	#print(pie[1,1])

	# use the recovered values to store the stopping probabilities
	for h=1:hours,sl=1:shiftl-1, c=1:num_types, a = 1:nwage

	    if h < hours 
	        nh = h+1
	    else 
	        nh = 1
	    end; 
	    ind = sl == shiftl - 1; 
	    nsl = min(sl + 1,shiftl-1);

	   	ANEV = 0.0
	    for ia = 1:nwage
	    	ANEV += pw[ia]*NEV[nh,nsl,c,ia]
	    end;
	    p_est[h,sl,c,a]=exp(1/σ)/(exp(1/σ)+exp((pie[h,a]*(1-breakm[h,sl,c,a])-λ_0[pick_o[h]]-λ_1*sl-λ_2*sl^2-conm[c,h]*f[pick_fine[h]]+(1-ind)*ANEV)/σ))
	
	end;

	# get starting probabilities
	for h=1:hours, c=1:num_types

		if h < hours 
	        nh = h+1
	    else 
	        nh = 1
	    end; 

		ANEV = 0.0
	    for ia = 1:nwage
	    	ANEV += pw[ia]*NEV[nh,1,c,ia]
	    end;

	    q_est[h,c]=exp((ANEV-lease[nh])/σ_s)/(exp((μ[work_time[nh]])/σ_s)+exp((ANEV-lease[nh])/σ_s))
	end;
	
	return p_est,q_est, NEV;	
end;


