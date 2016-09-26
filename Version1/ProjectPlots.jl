# Basic Plots

using DataFrames
using Distributions
using Gadfly

yearins = [ [x; findfirst(data[:,fipscodeloc], x); findlast(data[:,fipscodeloc], x ); unique( data[findfirst(data[:,fipscodeloc], x):findlast(data[:,fipscodeloc], x ) , yearloc]  ) ] for x in unique(data[:,fipscodeloc])  ];


monopoly = Array{Int64}(0)
duopoly = Array{Int64}(0)
triopoly = Array{Int64}(0)
tetrapoly = Array{Int64}(0)
nopoly = Array{Int64}(0)
all = Array{Int64}(0)

monopoly2005 = Array{Int64}(0)
duopoly2005 = Array{Int64}(0)
triopoly2005 = Array{Int64}(0)
tetrapoly2005 = Array{Int64}(0)
nopoly2005 = Array{Int64}(0)
all2005 = Array{Int64}(0)

for year in 1990:2012
  for el in yearins
    unqfids = [x for x in unique(data[ (data[:,fipscodeloc].==el[1])&(data[:,yearloc].==year) ,fidloc])]
    push!(all, size(unqfids,1))
    if year == 2005
      push!(all2005, size(unqfids,1))
    end
  end
end

# All Market Years:
plot(x=all, Geom.histogram, Guide.xlabel("Hospitals Per Market"), Guide.ylabel("Frequency"), Guide.title("Market Sizes"))
plot(x = all[all.>2], Geom.histogram, Guide.xlabel("Hospitals Per Market"), Guide.ylabel("Frequency"), Guide.title("Market Sizes Greater than Two Hospitals"))

# 2005 only

plot(x = all2005, Geom.histogram, Guide.xlabel("Hospitals Per Market"), Guide.ylabel("Frequency"), Guide.title("Market Sizes, 2005 Only"))

plot(x = all2005[all2005.>2], Geom.histogram, Guide.xlabel("Hospitals Per Market"), Guide.ylabel("Frequency"), Guide.title("Market Sizes Greater than 2, 2005 Only"))









###
