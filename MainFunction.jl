# Main function



function Mainfun(dataf::DataFrame, people::DataFrame, mkt_fips::Int64, year::Int64, modelparameters::Array{Float64, 2}, entryprobs::Array{Float64}, fids::Array{Int64}; nsims = 10)
			 # returns a dataframe unless converted
      numfids = size(fids)[1]
      outp = zeros(numfids, 183)
			#Arguments: Simulator(dataf::DataFrame, peoplesub::DataFrame, year::Int64, mkt_fips::Int64,  state_history::Array{Float64,2}, demandmodelparameters::Array{Float64, 2}; T = 100, start = 2)
      print("Equilibrium Simulation, ", mkt_fips, " ", year, " ", "\n")
                                      # Careful with this name
      states = Simulator(dataf, people, year, mkt_fips, modelparameters, entryprobs, T = 100, sim_start = 2)
      # Non-equilibrium Play -
			# Entrants in dataframe now tagged with negative ID's.  Remake to remove them:
			dataf = dataf[(dataf[:id].>= 0)&(!isna(dataf[:fipscode])), :];
		for f in 1:numfids
        pfid = fids[f]
        print("Perturbing Fid: ", pfid, "\n")
  			#Arguments: function PerturbSimulator(dataf::DataFrame, peoplesub::DataFrame, subname::ASCIIString, year::Int64, mkt_fips::Int64, demandmodelparameters::Array{Float64, 2}, pfid::Int64; disturb = 0.05, T = 100, sim_start = 2)
        perturbed_history = PerturbSimulator(dataf, peoplesub, year, mkt_fips, modelparameters, pfid, entryprobs, disturb = 0.01, T = 100, sim_start = 2)

  			# Here apply DynamicValue to the result of the simulations
  			# DynamicValue(state_history::Array, fac_fid::Float64; pat_types = 1, β = 0.95, T = 100, max_hosp = 25)
  			# output is in format: facility changes record, per-period visits record.
  			pfid_f = convert(Float64, pfid)
  			eq_change, eq_val  = DynamicValue(states, pfid_f; pat_types = 1, β = 0.95, T = 100, max_hosp = 25)
  			neq_change, neq_val = DynamicValue(perturbed_history, pfid_f; pat_types = 1, β = 0.95, T = 100, max_hosp = 25)
        outp[f,:] = [mkt_fips pfid_f year eq_val eq_change neq_val neq_change]
        # Abandon Entrants again.
        dataf = dataf[(dataf[:id].>= 0)&(!isna(dataf[:fipscode])), :]
    end
  return outp
end
