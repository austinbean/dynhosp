using DataFrames
using DataArrays
# There is a problem here: what probs should be returned when the hospital level ISN'T present?
function findprob2(dat1::DataFrame, Total::Int64, Lev1::Int64, Lev2::Int64, Lev3::Int64)
  if (Lev1 + Lev2 + Lev3 != Total)
    return "Level 1 + Level 2 + Level 3 must equal Total"
  else
      # Finds markets in the data frame matching some particular configuration of levels
      ind1 = dat1[:, :level1_hospitals0] .== Lev1
      ind2 = dat1[:, :level2solo_hospitals0] .== Lev2
      ind3 = dat1[:, :level3_hospitals0] .== Lev3
      ind4 = dat1[:, :total_hospital0] .== Total
      all_ind = (ind1.*ind2).*(ind4.*ind3)
      if sum(all_ind) == 0
        return "No such config"
      else
        # Picks out choices made by firms in that market configuration
        choices =sortrows( vcat( hcat( dat1[all_ind, :choicenum0], dat1[all_ind, :pr_ch_0], dat1[ all_ind , :act_int], dat1[ all_ind , :act_solo]), hcat(dat1[all_ind, :choicenum1], dat1[all_ind, :pr_ch_1], dat1[ all_ind , :act_int], dat1[ all_ind , :act_solo] ), hcat(dat1[all_ind, :choicenum2], dat1[all_ind, :pr_ch_2], dat1[ all_ind , :act_int], dat1[ all_ind , :act_solo]), hcat(dat1[all_ind, :choicenum3], dat1[all_ind, :pr_ch_3], dat1[ all_ind , :act_int], dat1[ all_ind , :act_solo])  ), by = x -> x[1])
        # "continue" choices differ: one for each level here.
        cont1 = choices[( (choices[:,1].==10).*(choices[:,3].==0).*(choices[:,4].==0) ), :  ] #note that the parenthese are crucial around the elements of the first index
        cont1 = cont1[findin(isna(cont1[:,2]), false), :]
        cont2 = choices[ ((choices[:,1].==10).*(choices[:,3].==0).*(choices[:,4].==1) ) ,:]
        cont2 = cont2[findin(isna(cont2[:,2]), false), :]
        cont3 = choices[ ((choices[:,1].==10).*(choices[:,3].==1).*(choices[:,4].==0)  ) ,:]
        cont3 = cont3[findin(isna(cont3[:,2]), false), :]
      #=
       1 "To 3 from 1" 2 "To 2 from 1" 3 "To 2 from 3" 4 "To 1 from 3" 5 "To 1 from 2" 6 "To 3 from 2" 7 "Enter at 1" 8 "Enter at 2" 9 "Enter at 3" 10 "Do Nothing" 11 "Exit"
      =#
          # Level 1
          p11  = maximum(min(cont1[:,2], 1))
          p12  = choices[findfirst(choices[:, 1], 2), 2]
          p13  = choices[findfirst(choices[:, 1], 1), 2]
          p1ex = choices[findfirst(choices[:, 1], 11), 2]
          # Level 2
          p21  = choices[findfirst(choices[:, 1], 5), 2]
          p22  = maximum(min(cont2[:,2], 1))
          p23  = choices[findfirst(choices[:, 1], 6), 2]
          p2ex = choices[findfirst(choices[:, 1], 11), 2]
          # Level 3
          p31  = choices[findfirst(choices[:, 1], 4), 2]
          p32  = choices[findfirst(choices[:, 1], 3), 2]
          p33  = maximum(min(cont3[:,2], 1))
          p3ex = choices[findfirst(choices[:, 1], 11), 2]
          return [p11, p12, p13, p1ex, p21, p22, p23, p2ex, p31, p32, p33, p3ex]
      end
    end
end
