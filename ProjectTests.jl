# ProjectTests.jl

using Base.Test

dyn = CounterObjects(50);


@testset "ProbUpdate Tests" begin 
  dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)] = nlrec( Dict(1=>1.0, 2=>1.0, 3=>0.0, 4=>0.0), zeros(2,4), Dict(99=>0) )
  ProbUpdate(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].aw, dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi)
  ix1 = findin(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi[1,:], 1)
  ix4 = findin(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi[1,:], 4)
  Base.Test.@test isapprox(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi[2,ix1], [0.36552928931500245])
  Base.Test.@test isapprox(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi[2,ix4], [0.13447071068499755])
  # change the values: 
  dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].aw[2] = 0;
  dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].aw[1] = 0;
  dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].aw[3] = 1;
  dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].aw[4] = 1;
  ProbUpdate(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].aw, dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi)
  Base.Test.@test isapprox(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi[2,ix1], [0.13447071068499755]) 
  Base.Test.@test isapprox(dyn.all[2].visited[(0,0,0,0,0,0,0,0,0,99)].psi[2,ix4], [0.36552928931500245]) 
end 


@testset "PolicyUpdate Tests:" begin
    pu1 = rand(2,10)
    PolicyUpdate(pu1)
    @test sum(pu1[2,:]) ≈ 1.0
end