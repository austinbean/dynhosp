using DataFrames

test_vals = zeros(Int, 1, 27 )
#=
 1 "To 3 from 1"
 2 "To 2 from 1"

 3 "To 2 from 3"
 4 "To 1 from 3"

 5 "To 1 from 2"
 6 "To 3 from 2"

 7 "Enter at 1"
 8 "Enter at 2"
 9 "Enter at 3"

 10 "Do Nothing"
 11 "Exit"
=#

count1 = 10000;
choices1 = [1, 2, 10, 11]; # for level 1
choices2 = [5, 6, 10, 11]; # for level 2
choices3 = [3, 4, 10, 11]; # for level 3
choices4 = [7, 8, 9, 10]; # for entrants
max_hosp = 22

# This one has a problem with the way level 2 and higher facilities are stored.
#=
Column Format:
[1 - fid, 2 - year, 3 - choicenum, 4 - ObstetricsLevel,
5 - NeoIntensive,
6 - TotalDeliveries, 7 - fipscode,
8 - SoloIntermediate,
9 - Choice, 10 - level3_entrants, 11 - level2_soloentrants,
12 - level3_exits, 13 - level2_soloexits,
14 - level3_hospitals,
15 - level2solo_hospitals,
16 - level1_hospitals,
17 - level3_potentials,
18 - level2solo_potentials,
19 - total_hospital,
20 - total_exits, 21 - nlevel1_entrants, 22 - nlevel2_entrants, 23 -nlevel3_entrants, 24 - active,
25 - act_int, 26 - act_solo, 27 - id]

=#

for i = 1:max_hosp # total hospitals
	println(i)
  for j = 0:i  # actual level 2's
    for k = 0:i # actual level 3's
			if j + k <= i  # hospital levels at 2, 3 must be less than total
				if (i-j-k) > 0
					for ch1 in choices1
						if ch1 == 1
							α = 1; β = 0; γ = -1; δ = -1; ϵ = 0; ζ = 0; χ = 1; ξ = 0;
						elseif ch1 == 2
							α = 0; β = 1; γ = -1; δ = 0; ϵ = -1; ζ = 0; χ = 0; ξ = 1;
						elseif ch1 == 10
							α = 0; β = 0; γ = 0; δ = 0; ϵ = 0; ζ = 0; χ = 0; ξ = 0;
						elseif ch1 == 11
							α = 0; β = 0; γ = -1; δ = -1; ϵ = -1; ζ = -1; χ = 0; ξ = 0;
						else
							break
						end
						count1 = count1 + 1           #column  5     8                 14            15         16                17                    18         19
						test_vals = vcat(test_vals, [0 0 ch1 0 χ 0 0 ξ 1 0 0 0 0 min((k+α),i) min((j+β), i) max((i-j-k + γ),0) max((i-k + δ),0) max((i-j + ϵ),0) (i +ζ) 0 0 0 0 ch1 0 0 count1])
							for x1 in setdiff(choices1, ch1)
								if x1 == 1
									α1 = 1; β1 = 0; γ1 = -1; δ1 = -1; ϵ1 = 0; ζ1 = 0; χ1 = 1; ξ1 = 0;
								elseif x1 == 2
									α1 = 0; β1 = 1; γ1 = -1; δ1 = 0; ϵ1 = -1; ζ1 = 0; χ1 = 0; ξ1 = 1;
								elseif x1 == 10
									α1 = 0; β1 = 0; γ1 = 0; δ1 = 0; ϵ1 = 0; ζ1 = 0; χ1 = 0; ξ1 = 0;
								elseif x1 == 11
									α1 = 0; β1 = 0; γ1 = -1; δ1 = -1; ϵ1 = -1; ζ1 = -1; χ1 = 0; ξ1 = 0;
								else
									break
								end
								test_vals = vcat(test_vals, [0 0 x1 0 χ1 0 0 ξ1 0 0 0 0 0 min((k+α1),i) min((j+β1), i) max((i-j-k + γ1),0) max((i-k + δ1),0) max((i-j + ϵ1),0) (i +ζ1) 0 0 0 0 ch1 0 0 count1])
							end
					end
				end
				if j > 0
					for ch2 in choices2
						if ch2 ==5
							α = 0; β = -1; γ = 1; δ = 0; ϵ = 1; ζ = 0; χ = 0; ξ = 0;
						elseif ch2 ==6
							α = 1; β = -1; γ = 0; δ = -1; ϵ = 1; ζ = 0; χ = 1; ξ = 0;
						elseif ch2 ==10
							α =0 ; β = 0; γ = 0; δ = 0 ; ϵ = 0; ζ = 0; χ = 0; ξ = 1;
						elseif ch2 ==11
							α = 0; β = -1; γ = 0; δ = -1; ϵ = 0; ζ = -1; χ = 0; ξ = 0;
						else
							break
						end
						count1 = count1 + 1                             #column        14            15         16                17                    18         19
						test_vals = vcat(test_vals, [0 0 ch2 0 χ 0 0 ξ 1 0 0 0 0 min((k+α),i) min((j+β), i) max((i-j-k + γ),0) max((i-k + δ),0) max((i-j + ϵ),0) (i +ζ) 0 0 0 0 ch2 0 1 count1])
						for x2 in setdiff(choices2, ch2)
							if x2== 5
								α1 = -1; β1 = -1; γ1 =1 ; δ1 = 0; ϵ1 = 1; ζ1 = 0; χ1 = 0; ξ1 = 0;
							elseif x2 == 6
								α1 = 1; β1 = -1; γ1 = 0; δ1 = -1; ϵ1 = 1; ζ1 = 0; χ1 = 1; ξ1 = 0;
							elseif x2 ==10
								α1 = 0; β1 = 0; γ1 = 0; δ1 = 0; ϵ1 = 0; ζ1 = 0; χ1 = 0; ξ1 = 1;
							elseif x2 ==11
								α1 = 0; β1 = -1; γ1 = 0; δ1 = -1; ϵ1 = 0; ζ1 = -1; χ1 = 0; ξ1 = 0;
							else
								break
							end
							test_vals = vcat(test_vals, [0 0 x2 0 χ1 0 0 ξ1 0 0 0 0 0 min(max((k+α1),0),i) min((j+β1), i) max((i-j-k + γ1),0) max((i-k + δ1),0) max((i-j + ϵ1),0) (i +ζ1) 0 0 0 0 ch2 0 1  count1])
						end
					end
				end
				if k > 0
					for ch3 in choices3
						if ch3 == 4 #note order here 4 and then 3
							α = -1; β = 0; γ = 1; δ =1 ; ϵ = 0; ζ = 0; χ = 0; ξ = 0;
						elseif ch3 == 3
							α = -1; β = 1; γ = 0; δ = 1; ϵ = -1; ζ = 0; χ = 0; ξ = 1;
						elseif ch3 == 10
							α = 0; β = 0; γ = 0; δ = 0; ϵ = 0; ζ = 0; χ = 1; ξ = 0;
						elseif ch3 == 11
							α = -1; β = 0; γ = 0; δ = 0; ϵ = -1; ζ = -1; χ = 0; ξ = 0;
						else
							break
						end
						count1 = count1 + 1                             #column        14            15         16                17                    18         19
						test_vals = vcat(test_vals, [0 0 ch3 0 χ 0 0 ξ 1 0 0 0 0 min((k+α),i) min((j+β), i) max((i-j-k + γ),0) max((i-k + δ),0) max((i-j + ϵ),0) (i +ζ) 0 0 0 0  ch3 1 0 count1])
						for x3 in setdiff(choices3, ch3)
							if x3 == 4 # note order here, 4 and then 3
								α1 = -1; β1 = 0; γ1 = 1; δ1 = 1; ϵ1 = 0; ζ1 = 0; χ1 = 0; ξ1 = 0;
							elseif x3 == 3
								α1 = -1; β1 = 1; γ1 = 0; δ1 = 1; ϵ1 = -1; ζ1 = 0; χ1 = 0; ξ1 = 1;
							elseif x3 == 10
								α1 = 0; β1 = 0; γ1 = 0; δ1 = 0; ϵ1 = 0; ζ1 = 0; χ1 = 1; ξ1 = 0;
							elseif x3 == 11
								α1 = -1; β1 = 0; γ1 = 0; δ1 = 0; ϵ1 = -1; ζ1 = -1; χ1 = 0; ξ1 = 0;
							else
								break
							end
							test_vals = vcat(test_vals, [0 0 x3 0 χ1 0 0 ξ1 0 0 0 0 0 min(max((k+α1),0),i) min((j+β1), i) max((i-j-k + γ1),0) max((i-k + δ1),0) max((i-j + ϵ1),0) (i +ζ1) 0 0 0 0 ch3 1 0 count1])
						end
					end
				end
				# for ch4 in choices4
				# 	count1 = count1 + 1
				# 	test_vals = vcat(test_vals, [0 0 ch4 0 0 0 0 0 1 0 0 0 0 k j (i-j-k) (i-k) (i-j) i 0 0 0 0 ch4 1 0 count1])
				# 	for x4 in setdiff(choices4, ch4)
				# 		test_vals = vcat(test_vals, [0 0 x4 0 0 0 0 0 0 0 0 0 0 k j (i-j-k) (i-k) (i-j) i 0 0 0 0 ch4 1 0 count1])
				# 	end
				# end
      end
    end
  end
end



df_simulated = DataFrame(test_vals[2:end, :])
d_columns = Dict{Any, Any}( :x1 => :fid, :x2 => :year, :x3 => :choicenum, :x4 => :ObstetricsLevel, :x5 => :NeoIntensive, :x6 => :TotalDeliveries, :x7 =>:fipscode, :x8 => :SoloIntermediate, :x9 => :choice, :x10 =>:level3_entrants, :x11 => :level2_soloentrants, :x12 => :level3_exits, :x13 => :level2_soloexits, :x14 => :level3_hospitals, :x15 => :level2solo_hospitals, :x16 => :level1_hospitals, :x17 => :level3_potentials, :x18 => :level2solo_potentials, :x19 => :total_hospital, :x20 => :total_exits, :x21 => :nlevel1_entrants, :x22 => :nlevel2_entrants, :x23 => :nlevel3_entrants, :x24 => :active, :x25 => :act_int, :x26 => :act_solo, :x27 => :id )
rename!(df_simulated, d_columns)

writetable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/simulatedchoices.csv", df_simulated)
