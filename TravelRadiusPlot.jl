# Plot travel radii to hospitals, on average.


using PyPlot

# Data
# This comes from the 2005 1st Quarter PUDF.dta file, running TX Merge Hospital Choices.do.  Section is (right now) lines 229-240

# [ 5.035 ; 7.85 ; 9.26; 10.17; 10.5; 11.3; 11.6; 12.06; 12.48; 12.7 ]


theta = [0:2pi/30:2pi;]
data = [ 5.035  7.85 9.26 10.17 10.5 11.3 11.6 12.06 12.48 12.7 ] # distances to i-th closest hospital, where i is the coordinate in that vector, so 9.26 miles is third closest
chosen = 10.91 # average distance to closest


fig = figure("dist_avail", figsize=(10, 10))
ax = axes(polar = "true")
title("Distances to Nearest Hospitals in Miles \n Measured from Zip Code Centroid")
plot(theta, chosen*ones(size(theta)), "r--", label = "Chosen")
plot( [pi/3], [10.91], "o")

plot(theta, data[1]*ones(size(theta)), "c-",  label = "Closest")
plot([pi/4], [data[1]], "o")

plot(theta, data[2]*ones(size(theta)), "b-")
plot([pi/5], [data[2]], "o")

plot(theta, data[3]*ones(size(theta)), "b-")
plot([3pi/5], [data[3]], "o")

plot(theta, data[5]*ones(size(theta)), "g-", label = "5th Farthest")
plot([3pi/4], [data[5]], "o")

plot(theta, data[6]*ones(size(theta)), "m-", label = "6th Farthest")
plot([2pi], [data[6]], "o")

plot(theta, data[10]*ones(size(theta)), "y-", label = "Farthest")
plot([pi/2], [data[10]], "o")


ax[:set_thetagrids]([])
#ax[:grid](false)
ax[:text](0, 0, "0", style="italic")
legend(loc = "upper left", bbox_to_anchor=(-0.1,1.1))
arrow(0, 0,pi/8.5 , 14)
annotate("Distance in\n Miles", xy=(0, 0), xytext=[3*pi/4, 1])
fig[:canvas][:draw]()

#=
for i in 1:10
  print( "theta, data[$i]*ones(size(theta)), ")
end
=#
