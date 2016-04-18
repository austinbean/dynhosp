# I don't know if I actually need this at all.  State History already
# includes  everything I want.


function DynamicValue(state_history::Array, beta::Float64, T::Int64)
  len, width = size(state_history)
  fac_section = state_history[:, width-4]
  agg_section = state_history[:, width-3:end-1]
  probs = state_history[:, end]

  mtotal = maximum(sum(agg_section[:,1:3 ], 2)) #maximum total hospitals.
  mlev1 = maximum(agg_section[:,1]) # maximum level 1
  mlev2 = maximum(agg_section[:,2]) # maximum level 2
  mlev3 = maximum(agg_section[:,3]) # maximum level 3

  fids = unique(convert(Array{Int}, state_history[1, 1:6:end-4] ))
  numfids = size(fids)[1]
  mkt_history = zeros(6*T ,numfids)

  for ind in size(fids)[1]
    index = findfirst(state_history[1,:], fid[ind])
    history_slice = state_history[:, index:index+5]




  end


end
