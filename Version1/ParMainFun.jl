# A parallelizable version of MainFunction.  Runs faster in a parallel for loop.
#=
- Should be run in a construction like @parallel (+) for i = 1:nsims... etc.
- This does one stupid thing: the fid, year and mkt_fips are appended to the beginning of the record.
- This means they are added together with the linkage operation (+) in the @parallel for
- Divide by nsims to recover the right values.

=#

function ParMainfun(dataf::Matrix, insured::Matrix, insuredmodelparameters::Array{Float64, 2}, medicaid::Matrix, medicaidmodelparameters::Array{Float64, 2}, mkt_fips::Int64, year::Int64, fids::Array{Int64};  idloc = 1, npers = 50)
			 # returns a dataframe unless converted
      numfids = maximum(size(fids))
      outp = zeros(numfids, ((2*78+12)*2)+3) # 78 for Private, 78 for Medicaid, 12 for level changes or exit, x 2 for equilibrium or not, + 3 for mkt, fid, year idents.
      for i = 1:numfids
        outp[i,1] = mkt_fips
        outp[i,2] = fids[i]
        outp[i,3] = year
      end
      states =
			dataf = dataf[(dataf[:,idloc].>= 0), :];
      states = ProjectModule.Simulator(dataf, insured, insuredmodelparameters, medicaid, medicaidmodelparameters, year, mkt_fips; T = npers)
  		for f in 1:numfids
          pfid = fids[f]
          pfid_f = convert(Float64, pfid)
          # Output of DynamicValue is Medicaid Patients, level changes, Private Patients
          neq_med, neq_change, neq_val = ProjectModule.DynamicValue(ProjectModule.PerturbSimulator(dataf, insured, insuredmodelparameters, medicaid, medicaidmodelparameters, year, mkt_fips, pfid; disturb = 0.01, T = npers), pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
          eq_med, eq_change, eq_val = ProjectModule.DynamicValue(states, pfid_f; pat_types = 1, β = 0.95, max_hosp = 25)
          # Abandon Entrants again.
          dataf = dataf[(dataf[:,idloc].>= 0), :];
          outp[f,4:end] += [eq_val eq_med eq_change neq_val neq_med neq_change]
      end
  return outp
end
