# A parallelizable version of MainFunction.  Runs faster in a parallel for loop.
#=
- Should be run in a construction like @parallel (+) for i = 1:nsims... etc.
- This does one stupid thing: the fid, year and mkt_fips are appended to the beginning of the record.
- This means they are added together with the linkage operation (+) in the @parallel for
- Divide by nsims to recover the right values.

=#

function ParMainfun(dataf::Matrix, people::Matrix, mkt_fips::Int64, year::Int64, modelparameters::Array{Float64, 2}, fids::Array{Int64};  idloc = 1, npers = 50)
			 # returns a dataframe unless converted
      numfids = maximum(size(fids))
      outp = zeros(numfids, 183)
      for i = 1:numfids
        outp[i,1] = mkt_fips
        outp[i,2] = fids[i]
        outp[i,3] = year
      end
      states =
			dataf = dataf[(dataf[:,idloc].>= 0), :];
      states = ProjectModule.Simulator(dataf, people, year, mkt_fips, modelparameters; T = npers)
  		for f in 1:numfids
          pfid = fids[f]
          pfid_f = convert(Float64, pfid)
          # Note that it is important that all user-defined functions in this function say ProjectModule.(fn) - these are imported by the module but not brought into scope.
          # This probably won't be necessary if this is a function in its own file and imported directly via the module.  (The other two functions would surely cause problems if that were necessary.)
          neq_change, neq_val = ProjectModule.DynamicValue(ProjectModule.PerturbSimulator(dataf, people, year, mkt_fips, modelparameters, pfid; disturb = 0.01, T = npers), pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
          eq_change, eq_val = ProjectModule.DynamicValue(states, pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
          # Abandon Entrants again.
          dataf = dataf[(dataf[:,idloc].>= 0), :];
          outp[f,4:end] += [eq_val eq_change neq_val neq_change]
      end
  return outp
end
