# This tests the variance of demand against the variance from the data.
#=
Structure -
- Take demand Model
- Simulate nsims years
- There should be variation driven by the random shock
- Compare this to the variation in total deliveries reported.
- Ponder what to do next.

=#

include("/Users/austinbean/Desktop/dynhosp/Reboot.jl")
using Gadfly
using DataFrames

empirical_sds = DataFrames.readtable("/Users/austinbean/Desktop/dynhosp/birthvariances.csv", header = true);

A = unique(pinsured[:,2])
const dims1 = length(A)
const nsims = 50
private_vals = zeros(dims1, nsims+1)
medicaid_vals = zeros(dims1, nsims+1)
totals = zeros(dims1, nsims+1)
set1 = collect(unique(pinsured[:,2]))
for el in 1:length(set1)
  totals[el,1] = set1[el]
  private_vals[el,1] = set1[el]
  medicaid_vals[el,1] = set1[el]
end

for i = 1:nsims
  priv = countmap(DemandModel(pinsured, privatedemandmodelparameters, Array{Float64,2}()))
  for j in keys(priv)
    if j > 0
      index = findfirst(private_vals[:,1], j)
      if index > 0
    #    println(priv[j])
        private_vals[index, i+1] += priv[j]
      end
    end
  end

  med = countmap(DemandModel(pmedicaid, medicaiddemandmodelparameters, Array{Float64,2}()))
  for j in keys(med)
      index = findfirst(medicaid_vals[:,1], j )
      if index > 0
    #    println(med[j])
        medicaid_vals[index,i+1] += med[j]
      end
    end
end

for i = 1:size(private_vals,1)
  for j = 2:size(private_vals,2)
    totals[i,j] += private_vals[i,j]
    totals[i,j] += medicaid_vals[i,j]
  end
end

sds = zeros(dims1,2)

for i=1:dims1
  sds[i,1] = totals[i,1]
  sds[i,2] = std( totals[i, 2:end])
end

plot(layer(x=sds[ :,2], Geom.histogram, Theme(default_color=colorant"blue")),
     layer(x = empirical_sds[ empirical_sds[:sdd].<1000, :sdd], Geom.histogram, Theme(default_color=colorant"red")),
     Guide.xlabel("Variance"), Guide.ylabel("Frequency"), Guide.manual_color_key("Models", ["Choice Model", "Data"], ["blue",  "red"]))

# multiply the variance in sds by 10 and the plot looks much better.





###
