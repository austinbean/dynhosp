

function DynamicValue(state_history::Array, fac_fid::Float64; pat_types = 1, β = 0.95, T = 100, max_hosp = 25)
  len, width = size(state_history)
  index = findfirst(state_history[1,:], fac_fid) # where is the perturbed facility
  fac_section = state_history[:, index:index+6] # take the section from the actual history corresponding
  agg_section = state_history[:, width-3:end]
  history = hcat(fac_section, agg_section)

#= Output format:
 Row 1 : [0, 1, 2, ..., 25] - for level 1
 Row 2 : [0, 1, 2, ..., 25] - for level 2
 Row 3 : [0, 1, 2, ..., 25] - for level 3

- Each row element will be: Β^t * Demand * probability of visit
- The choice of row depends on the level of the hospital's facility
- The specific element of each row depends on how many facilities at the level there
  are
- Think ahead: if we do this separately for Medicaid patients, I'll need 6 rows since the per-patient
  revenues are different
- Entry and exit to be recorded in a separate manner.
=#

  outp = zeros(3*pat_types, max_hosp + 1) # visits to max hospital or 0 at each of 3 levels

# The last 4 columns of "history" or "state_history" are: level 1, level 2, level 3, probability
# The seven elements of the records are: fid, act_solo, act_int, choice prob, action, demand and perturbation indicator

  for row in 1:(T+1):
    if (history[row,2], history[row,3]) == (0,0)
      levelcount = history[row, end-3]
      outp[3, levelcount+1] += (β^row)*history[row, 6]*history[row,end]
    elseif (history[row,2], history[row,3]) == (1,0)
      levelcount = history[row, end-2]
      outp[2, levelcount+1] += (β^row)*history[row, 6]*history[row,end]
    elseif (history[row,2], history[row,3]) == (0,1)
      levelcount = history[row, end-1]
      outp[1, levelcount+1] += (β^row)*history[row, 6]*history[row,end]
    end
  end

return outp

end
