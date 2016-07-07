# DynamicValue.jl
# Takes the output of the simulation and stores that in "basis function" form for BBL

#= First Output Format - state visits, demand, probability:
 Row 1 : [0, 1, 2, ..., 25] - for level 1
 Row 2 : [0, 1, 2, ..., 25] - for level 2
 Row 3 : [0, 1, 2, ..., 25] - for level 3

- Each row element will be: Β^t * Demand * probability of visit
- The choice of row depends on the level of the hospital's facility
- The specific element of each row depends on how many facilities at the level there
  are
- Think ahead: if we do this separately for Medicaid patients, I'll need 6 rows since the per-patient
  revenues are different

Second output type - record facility changes:
- Return a vector [(12), (13), (1EX), (21), (23), (2EX), (31), (32), (3EX), Enter1, Enter2, Enter3]
- (1i) counts number of transitions from type 1 to i = 2, 3, EX
- (2i) counts number of transitions from type 2 to i = 1, 3, EX
- (3i) counts number of transitions from type 3 to i = 1, 2, EX
- Enter1 records entering at 1

- Thinking ahead: there may be differences depending on whether NFP or FP
- Thinking ahead: the 3 to 2 transition can also be a different sign (maybe hospital earns money selling off capital)
  than the 1 to 2 transition, which is surely negative.

  α₂ = 0.07; α₃ = 0.13; pat_types = 1; β = 0.95; max_hosp = 25
=#


function DynamicValue(state_history::Array, fac_fid::Float64; α₂ = 0.07, α₃ = 0.13, pat_types = 1, β = 0.95, max_hosp = 25)
  T, width = size(state_history)
  index = findfirst(state_history[1,:], fac_fid) # where is the perturbed facility
  #fac_section = state_history[:, index:index+6] # take the section from the actual history corresponding
  #agg_section = state_history[:, width-3:end]
  #history = hcat(fac_section, agg_section)
  history = hcat(state_history[:, index:index+6], state_history[:, width-3:end])
  if maximum(sum(state_history[:, end-3:end-1],2)) > 25
    println("Maximum hospitals observed exceeds max_hosp")
    return "Size Warning"
  end
  outp = zeros(3*pat_types, max_hosp + 1) # visits to max hospital or 0 at each of 3 levels
  outp2 = zeros(1, 12) # records

# The last 4 columns of "history" or "state_history" are: level 1, level 2, level 3, probability
# The seven elements of "history" (per fid in state_history) are:
# fid, act_solo, act_int, choice prob, action, XXXXX, demand and perturbation indicator

# Write out values for first row separately.

  if (history[1,2], history[1,3]) == (0,0) # hosp starts at level 1
    levelcount = convert(Int64, sum(history[1, end-3:end-1])) # counts the number of hospitals total
    outp[1,levelcount+1] += (β^0)*history[1, 6]*history[end,end]
  elseif (history[1,2], history[1,3]) == (1,0)
    levelcount = convert(Int64, sum(history[1, end-3:end-1])) # count of all facilities
    levelcount2 = convert(Int64, history[1,end-2]) # count of level 2 facilities
    outp[1, levelcount+1] += (β^0)*history[1, 6]*history[end,end] # regular births
    outp[2, levelcount2+1] += (α₂)*(β^0)*history[1, 6]*history[end,end] # expected fraction to NICU 2
  elseif (history[1,2], history[1,3]) == (0,1)
    levelcount = convert(Int64, sum(history[1, end-3:end-1]))
    levelcount3 = convert(Int64, history[1,end-1]) # count of level 3's
    outp[1, levelcount+1] += (β^0)*history[1, 6]*history[end,end]
    outp[3, levelcount3+1] += (α₃)*(β^0)*history[1, 6]*history[end,end] # expected fraction to NICU 3
  elseif (history[1,2], history[1,3]) == (-999,-999)
    # This is a firm which exits as a first action
    # println( "Firm is exiting in first period - ?") # need to think about what to do about this possibility.
  elseif (history[1,2], history[1,3]) == (999,999)
    # do nothing with this firm - it enters later.
  end
  for row in 2:T
    if (history[row,2], history[row,3]) == (0,0)
      levelcount = convert(Int64, sum(history[row, end-3:end-1])) # revenue depends on total hospitals
      outp[1, levelcount+1] += (β^(row-1))*history[row, 6]*history[end,end] # this is discount^t * demand * probability (aggregate)
      if (history[row-1,2], history[row-1,3]) == (0,0)
        # do nothing
      elseif (history[row-1,2], history[row-1,3]) == (1,0)
        # here you downgraded 2 to 1
        outp2[1, 4] += 1*β^(row-1)*history[end,end]
      elseif (history[row-1,2], history[row-1,3]) == (0,1)
        # here you downgraded 3 to 1
        outp2[1, 7] += 1*β^(row-1)*history[end,end]
      elseif ((history[row-1,2], history[row-1,3]) == (999,999))
        # entered previous period
        outp2[1, 10] += 1*history[end,end]
      end
    elseif (history[row,2], history[row,3]) == (1,0)
      # Here is an error - end-2 on the next line.
      levelcount = convert(Int64, sum(history[row, end-3:end-1]))
      levelcount2 = convert(Int64, history[row,end-2]) # number of level 2's
      outp[2, levelcount+1] += (β^row)*history[row, 6]*history[end,end]
      outp[2, levelcount2+1] += (α₂)*(β^row)*history[row, 6]*history[end,end] # expected rev from lev 2 admissions
      if (history[row-1,2], history[row-1,3]) == (0,0)
        # Upgraded 1 to 2
        outp2[1,1] += 1*β^(row-1)*history[end,end]
      elseif (history[row-1,2], history[row-1,3]) == (1,0)
        # do nothing
      elseif (history[row-1,2], history[row-1,3]) == (0,1)
        # here you downgraded 3 to 2
        outp2[1,8] += 1*β^(row-1)*history[end,end]
      elseif ((history[row-1,2], history[row-1,3]) == (999,999))
        # entered previous period
        outp2[1,end-1] += 1*history[end,end]
      end
    elseif (history[row,2], history[row,3]) == (0,1)
      levelcount = convert(Int64, sum(history[row, end-3:end-1]))
      levelcount3 = convert(Int64, history[row,end-1])
      outp[3, levelcount+1] += (β^(row-1))*history[row, 6]*history[end,end]
      outp[3, levelcount3+1] += (α₃)*(β^(row-1))*history[row, 6]*history[end,end]
      if (history[row-1,2], history[row-1,3]) == (0,0)
        # upgraded 1 to 3
        outp2[1,2] += 1*(β^(row-1))*history[end,end]
      elseif (history[row-1,2], history[row-1,3]) == (1,0)
        #upgraded 2 to 3
        outp2[1,5] += 1*(β^(row-1))*history[end,end]
      elseif (history[row-1,2], history[row-1,3]) == (0,1)
        # do nothing
      elseif ((history[row-1,2], history[row-1,3]) == (999,999))
        # entered previous period
        outp2[1,end] += 1*history[end,end]
      end
    elseif (history[row,2], history[row,3]) == (-999,-999)
      # here the firm has exited - need to record this the first time only
      # do this separately.
      if (history[row-1,2], history[row-1,3]) == (0,0)
        # Exited at 1
        outp2[1,3] += 1*β^(row-1)*history[end,end]
      elseif (history[row-1,2], history[row-1,3]) == (1,0)
        # Exited at 2
        outp2[1,6] += 1*β^(row-1)*history[end,end]
      elseif (history[row-1,2], history[row-1,3]) == (0,1)
        # Exited at 3
        outp2[1,9] += 1*β^(row-1)*history[end,end]
      elseif ((history[row-1,2], history[row-1,3]) == (999,999))
        # entered previous period
         outp2[1,3] += 0 # this is a guy who enters and immediately exits.  Not a case worth worrying about, I think.
      end
    else
      println("at row ", row, " solo ", history[row,2], " intensive ", history[row,3] )
      return "Bad firm state"
    end
  end
  outp = [ outp[1,1:end] outp[2,1:end] outp[3,1:end]] #rearranges the matrix.
return outp2, outp
end



#=
# Testing:

# This should work and generate 10*probability
state_history2 = [111 1 0 1/2 10 10 0 1 1 1 1/2 ; 111 0 1 1/3 10 10 0 1 1 1 1/3 ; 111 0 0 1/4 10 10 0 1 1 1 1/4 ]
invest, hist = DynamicValue(state_history2, 111.0, β= 1)'

# This should give us 10*probability in different columns
state_history3 = [111 1 0 1/2 10 10 0 1 1 8 1/2 ; 111 0 1 1/3 10 10 0 1 9 1 1/3 ; 111 -999 -999 1/4 10 10 0 10 1 1 1/4 ]
invest3, hist3 = DynamicValue(state_history3, 111.0, β= 1)'

# This should give an error::
state_history = [111 1 0 1/2 10 10 0 1 1 25 1/2 ; 111 0 1 1/3 10 10 0 1 9 1 1/3 ; 111 0 0 1/4 10 10 0 10 1 1 1/4 ]
invest, hist = DynamicValue(state_history, 111.0, β= 1, T = 2)'

# Should record several entry and status changes.
state_history = [111 999 999 1/2 10 10 0 1 1 1 1/2 ; 111 0 0 1/2 10 10 0 1 1 1 1/2 ; 111 0 1 1/3 10 10 0 1 9 1 1/3 ; 111 0 0 1/4 10 10 0 10 1 1 1/4; 111 1 0 1/4 10 10 0 10 1 1 1/4; 111 0 1 1/4 10 10 0 10 1 1 1/4; 111 -999 -999 1/4 10 10 0 10 1 1 1/4 ]
invest, hist = DynamicValue(state_history, 111.0, β=1, T = 6)
=#
