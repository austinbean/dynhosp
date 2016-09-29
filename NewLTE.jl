# Redoing the LTE.

using DataFrames
using Distributions
using StatsBase
using Gadfly  # probably should be switched out for Plots


dat = readtable("/Users/austinbean/Google Drive/Simulation Results/combinedresults.csv");
