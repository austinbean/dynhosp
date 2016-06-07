# Main function



function Mainfun(dataf::Matrix, people::Matrix, mkt_fips::Int64, year::Int64, modelparameters::Array{Float64, 2}, fids::Array{Int64}; nsims = 500)
			 # returns a dataframe unless converted
      numfids = size(fids)[1]
      outp = zeros(numfids, 183)
      for j in 1:nsims
  			#Arguments: Simulator(data::Matrix, peoplesub::Matrix, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2};  T = 100, sim_start = 2, fields = 7, neighbors_start = 108, entrants = [0, 1, 2, 3], entryprobs = [0.9895, 0.008, 0.0005, 0.002])
    #    print("Equilibrium Simulation, ", mkt_fips, " ", year, " ", "\n")
        states = Simulator(dataf, people, year, mkt_fips, modelparameters; T = 2)
        # Non-equilibrium Play -
  			# Entrants in dataframe now tagged with negative ID's.  Remake to remove them:
  			dataf = dataf[(dataf[:,idloc].>= 0), :];
    		for f in 1:numfids
            pfid = fids[f]
        #    print("Perturbing Fid: ", pfid, "\n")
      			#Arguments: function PerturbSimulator(data::Matrix, peoplesub::Matrix, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2}, pfid::Int64; entryprobs = [0.9895, 0.008, 0.0005, 0.002], entrants = [0, 1, 2, 3], disturb = 0.05,  T = 100, sim_start = 2, fields = 7, neighbors_start = 108)
# Here is a place to parallelize across the different fids
            perturbed_history = PerturbSimulator(dataf, people, year, mkt_fips, modelparameters, pfid; disturb = 0.01, T = 2)

      			# Here apply DynamicValue to the result of the simulations
      			# DynamicValue(state_history::Array, fac_fid::Float64; α₂ = 0.07, α₃ = 0.13, pat_types = 1, β = 0.95, max_hosp = 25)
      			# output is in format: facility changes record, per-period visits record.
      			pfid_f = convert(Float64, pfid)
      			eq_change, eq_val  = DynamicValue(states, pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
      			neq_change, neq_val = DynamicValue(perturbed_history, pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
            if j == 1
              outp[f,:] = [mkt_fips pfid_f year eq_val eq_change neq_val neq_change]
            else
              outp[f,4:end] = outp[f,4:end] + [eq_val eq_change neq_val neq_change]
            end
            # Abandon Entrants again.
            dataf = dataf[(dataf[:,idloc].>= 0), :];
        end
    end
  return outp
end
