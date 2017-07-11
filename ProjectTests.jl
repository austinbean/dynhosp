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


@testset "ConvTest tests" begin
    ab = Dict( 1 => true, 2 => true, 3 =>true)  # three firms converged.  Return false.  
    @test !ConvTest(ab) 

    ab2 = Dict( 1 => true, 2 => true, 3 =>false) # two firms converged, one not.  Return true.
    @test ConvTest(ab2) 

    ab3 = Dict(1 => false, 2=> false, 3 => false) # START - no one has converged.  Return true.
    @test ConvTest(ab3) 
end


@testset "FindComps Test" begin 
    ns = Array{Int64,1}()
    FindComps(dyn, ns , dyn.all[11]) 
    @test ns == [195, 196, 197, 198]
    FindComps(dyn, ns , dyn.all[11]) 
    @test ns == [195, 196, 197, 198] # call twice to make sure it doesn't append the same numbers.  
end 

@testset "MakeStateBlock Test" begin 

    @test MakeStateBlock([10]) == [(10, 1);(10, 2);(10, 3);(10, 999)]

end 


@testset "NFids test" begin 
    nf = Array{Int64,1}()
    NFids(dyn, nf, dyn.all[1])
    @test nf == [1391330]
    NFids(dyn, nf, dyn.all[1])
    @test nf == [1391330] # make sure it does not reappend.  
    ns = Array{Int64,1}()
    NFids(dyn, ns, dyn.all[11], dyn.all[13])
    @test ns == [3396057, 3390720, 3396327, 3396189,616303, 610460, 616318, 611790]
end 


@testset "ResultsOutVariant" begin 

  Texas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
  Nexas = CreateEmpty(ProjectModule.fips, ProjectModule.alldists, 50);
  Restore(Texas) # assign all to 0 or 1
  Restore(Nexas)

##  Test ##
  Texas.mkts[48453].config[1].levelhistory[2] = 2
  Texas.mkts[48453].config[1].levelhistory[3] = 3 
  Texas.mkts[48453].config[1].wtphist.w385[1] = 1.0
  Texas.mkts[48453].config[1].wtphist.w386[1] = 1.0
  Texas.mkts[48453].config[1].wtphist.w387[1] = 1.0
  Texas.mkts[48453].config[1].wtphist.w388[1] = 1.0
  Texas.mkts[48453].config[1].wtphist.w389[1] = 1.0
  Texas.mkts[48453].config[1].wtphist.w390[1] = 1.0
  Texas.mkts[48453].config[1].wtphist.w391[1] = 1.0
  Texas.mkts[48453].config[1].pdemandhist.demand385[1] = 1
  Texas.mkts[48453].config[1].pdemandhist.demand386[1] = 1
  Texas.mkts[48453].config[1].pdemandhist.demand387[1] = 1
  Texas.mkts[48453].config[1].pdemandhist.demand388[1] = 1
  Texas.mkts[48453].config[1].pdemandhist.demand389[1] = 1
  Texas.mkts[48453].config[1].pdemandhist.demand390[1] = 1
  Texas.mkts[48453].config[1].pdemandhist.demand391[1] = 1
  Texas.mkts[48453].config[1].mdemandhist.demand385[1] = 1
  Texas.mkts[48453].config[1].mdemandhist.demand386[1] = 1
  Texas.mkts[48453].config[1].mdemandhist.demand387[1] = 1
  Texas.mkts[48453].config[1].mdemandhist.demand388[1] = 1
  Texas.mkts[48453].config[1].mdemandhist.demand389[1] = 1
  Texas.mkts[48453].config[1].mdemandhist.demand390[1] = 1
  Texas.mkts[48453].config[1].mdemandhist.demand391[1] = 1
    # noneq 
  Nexas.mkts[48453].config[1].levelhistory[2] = 2
  Nexas.mkts[48453].config[1].levelhistory[3] = 3 
  Nexas.mkts[48453].config[1].levelhistory[4] = 1
  Nexas.mkts[48453].config[1].levelhistory[5] = 3
  Nexas.mkts[48453].config[1].levelhistory[6] = 2
  Nexas.mkts[48453].config[1].levelhistory[7] = 1
  Nexas.mkts[48453].config[1].levelhistory[8] = -999
  Nexas.mkts[48453].config[1].wtphist.w385[1] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w386[1] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w387[1] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w388[1] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w389[1] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w390[1] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w391[1] = 1.0
    # for level 2 
  Nexas.mkts[48453].config[1].wtphist.w385[2] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w386[2] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w387[2] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w388[2] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w389[2] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w390[2] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w391[2] = 1.0
    # for level 3 
  Nexas.mkts[48453].config[1].wtphist.w385[3] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w386[3] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w387[3] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w388[3] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w389[3] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w390[3] = 1.0
  Nexas.mkts[48453].config[1].wtphist.w391[3] = 1.0
    # demand level 1 
  Nexas.mkts[48453].config[1].pdemandhist.demand385[1] = 1
  Nexas.mkts[48453].config[1].pdemandhist.demand386[1] = 1
  Nexas.mkts[48453].config[1].pdemandhist.demand387[1] = 1
  Nexas.mkts[48453].config[1].pdemandhist.demand388[1] = 1
  Nexas.mkts[48453].config[1].pdemandhist.demand389[1] = 1
  Nexas.mkts[48453].config[1].pdemandhist.demand390[1] = 1
  Nexas.mkts[48453].config[1].pdemandhist.demand391[1] = 1
  Nexas.mkts[48453].config[1].mdemandhist.demand385[1] = 1
  Nexas.mkts[48453].config[1].mdemandhist.demand386[1] = 1
  Nexas.mkts[48453].config[1].mdemandhist.demand387[1] = 1
  Nexas.mkts[48453].config[1].mdemandhist.demand388[1] = 1
  Nexas.mkts[48453].config[1].mdemandhist.demand389[1] = 1
  Nexas.mkts[48453].config[1].mdemandhist.demand390[1] = 1
  Nexas.mkts[48453].config[1].mdemandhist.demand391[1] = 1
    

  ans1 = ResultsOutVariant(Texas, Texas; T = 0) # no discounting , beta = 0 must be the same.
  indx = findin(ans1[:,1], 4536253) # 261
  isapprox(ans1[indx,2:23],ans1[indx,24:end])   # one test 
  ans2 = ResultsOutVariant(Texas, Texas; beta = 0.0)
  ans3 = ResultsOutVariant(Nexas,Nexas;T=0)
  
  # Actual Tests. 
  Base.Test.@test isapprox(ans3[261,2],7.0)         # is the sum of WTP and demand, all of which are ones, so should be 7 
  Base.Test.@test isapprox(ans3[261,3],0.0)         # is the sum of WTP and demand, all of which are ones, so should be 7
  Base.Test.@test isapprox(ans3[261,4],0.0)         # is the sum of WTP and demand, all of which are ones, so should be 7
  Base.Test.@test isapprox(ans3[261,5], 30.86)      # should be 2*sum weights.

  Base.Test.@test isapprox(ans3[261,8], 12038.83)   # Should be the DRG Amount for 385
  Base.Test.@test isapprox(ans3[261,9], 66143.19)   # Should be the DRG Amount for 386
  Base.Test.@test isapprox(ans3[261,10], 19799.52)  # Should be the DRG Amount for 387
  Base.Test.@test isapprox(ans3[261,11], 4044.67)   # Should be the DRG Amount for 388
  Base.Test.@test isapprox(ans3[261,12], 6242.39)   # Should be the DRG Amount for 389
  Base.Test.@test isapprox(ans3[261,13], 1329.98)   # Should be the DRG Amount for 390
  Base.Test.@test isapprox(ans3[261,14], 412.04)    # Should be the DRG Amount for 391
  Base.Test.@test isapprox(ans3[261,15], 1.0)       # Transition 1 to 2
  Base.Test.@test isapprox(ans3[261,16], 1.0)       # T 1 to 3
  Base.Test.@test isapprox(ans3[261,17], 1.0)       # T 1 to Ex 
  Base.Test.@test isapprox(ans3[261,18], 1.0)       # TODO 
  Base.Test.@test isapprox(ans3[261,19], 1.0)       # Transition 2 to 3


  # demand level 2
  Nexas.mkts[48453].config[1].pdemandhist.demand385[2] = 2
  Nexas.mkts[48453].config[1].pdemandhist.demand386[2] = 2
  Nexas.mkts[48453].config[1].pdemandhist.demand387[2] = 2
  Nexas.mkts[48453].config[1].pdemandhist.demand388[2] = 2
  Nexas.mkts[48453].config[1].pdemandhist.demand389[2] = 2
  Nexas.mkts[48453].config[1].pdemandhist.demand390[2] = 2
  Nexas.mkts[48453].config[1].pdemandhist.demand391[2] = 2
  Nexas.mkts[48453].config[1].mdemandhist.demand385[2] = 2
  Nexas.mkts[48453].config[1].mdemandhist.demand386[2] = 2
  Nexas.mkts[48453].config[1].mdemandhist.demand387[2] = 2
  Nexas.mkts[48453].config[1].mdemandhist.demand388[2] = 2
  Nexas.mkts[48453].config[1].mdemandhist.demand389[2] = 2
  Nexas.mkts[48453].config[1].mdemandhist.demand390[2] = 2
  Nexas.mkts[48453].config[1].mdemandhist.demand391[2] = 2
    # demand level 3 
  Nexas.mkts[48453].config[1].pdemandhist.demand385[3] = 3
  Nexas.mkts[48453].config[1].pdemandhist.demand386[3] = 3
  Nexas.mkts[48453].config[1].pdemandhist.demand387[3] = 3
  Nexas.mkts[48453].config[1].pdemandhist.demand388[3] = 3
  Nexas.mkts[48453].config[1].pdemandhist.demand389[3] = 3
  Nexas.mkts[48453].config[1].pdemandhist.demand390[3] = 3
  Nexas.mkts[48453].config[1].pdemandhist.demand391[3] = 3
  Nexas.mkts[48453].config[1].mdemandhist.demand385[3] = 3
  Nexas.mkts[48453].config[1].mdemandhist.demand386[3] = 3
  Nexas.mkts[48453].config[1].mdemandhist.demand387[3] = 3
  Nexas.mkts[48453].config[1].mdemandhist.demand388[3] = 3
  Nexas.mkts[48453].config[1].mdemandhist.demand389[3] = 3
  Nexas.mkts[48453].config[1].mdemandhist.demand390[3] = 3
  Nexas.mkts[48453].config[1].mdemandhist.demand391[3] = 3

  ans4 = ResultsOutVariant(Nexas,Nexas;T=0)

  Base.Test.@test isapprox(ans4[261,6], 61.72)      # should be 4*sum weights.
  Base.Test.@test isapprox(ans4[261,7], 92.58)      # should be 6*sum weights. 

  Texas = 0; Nexas = 0 # dump.  
  end # of begin block for test.  



