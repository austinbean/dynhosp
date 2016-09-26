# Model attempt

#=
These require solvers - some are free and some are not.
Gurobi can be downloaded for free.
CPLEX can be obtained from IBM, probably for free
MOSEK has a free academic version, renewable yearly
IPOPT seems to be free
KNITRO - 300 variables, 300 constraints 6 month license, maybe full version for a
one month license.

=#

using JuMP

gmm = Model()
@defVar(gmm, 0 <= x <= 10)
@defVar(gmm, 5 <= β <= 20)

@setObjective(gmm, :Max, 5x + 13β - (x + β)^2 )
@addConstraint(gmm, 1x + 1β <= 15.0)

print(gmm)

status = solve(gmm)

Pkg.add("Gurobi")
Pkg.add("CPLEX")
Pkg.add("Mosek")
Pkg.add("Ipopt")

Pkg.build("Gurobi")
