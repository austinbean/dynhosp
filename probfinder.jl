# Probfinder2


using DataFrames
using DataArrays
type ProbException <: Exception end


function findprob(dat1::DataFrame, Total::Int64, Lev1::Int64, Lev2::Int64, Lev3::Int64)
  if (Lev1 + Lev2 + Lev3 != Total)
    return "Level 1 + Level 2 + Level 3 must equal Total"
  else
      # Level 1's should be identified by the choices they have available, not facilities.
      #= Another trick - a market in which someone updated to Lev 3 has that value above 0.  If someone picked it *in that period* then Lev3 above will be greater than 0, but there will be NO firm with
      that set of choices.  So the set of indices below will be empty, which is correct, but the function's going to look for it.
      Solution:
      - there will be times when it definitely matters.
      - Look for the probs based on whether sum ind1 > 0
      - If it's 1 or greater, we know what to do.
      - When these things are empty, then return a prob error.  Then look in the simulated data.
      - That should do it.

         =#
        ind1 = (( (dat1[:,:choicenum0].==1) | (dat1[:,:choicenum0].==2) | (dat1[:,:choicenum1].==1) | (dat1[:,:choicenum1].==2) | (dat1[:,:choicenum2].==1) | (dat1[:,:choicenum2].==2) | (dat1[:,:choicenum3].==1) | (dat1[:,:choicenum3].==2)) & (dat1[:, :level1_hospitals0] .== Lev1 ) & (dat1[:, :level2solo_hospitals0] .== Lev2) & (dat1[:, :level3_hospitals0] .== Lev3))
        ind2 = (( (dat1[:,:choicenum0].==5) | (dat1[:,:choicenum0].==6) | (dat1[:,:choicenum1].==5) | (dat1[:,:choicenum1].==6) | (dat1[:,:choicenum2].==5) | (dat1[:,:choicenum2].==6) | (dat1[:,:choicenum3].==5) | (dat1[:,:choicenum3].==6)) & (dat1[:, :level1_hospitals0] .== Lev1 ) & (dat1[:, :level2solo_hospitals0] .== Lev2) & (dat1[:, :level3_hospitals0] .== Lev3))
        ind3 = (( (dat1[:,:choicenum0].==3) | (dat1[:,:choicenum0].==4) | (dat1[:,:choicenum1].==3) | (dat1[:,:choicenum0].==4) | (dat1[:,:choicenum2].==3) | (dat1[:,:choicenum0].==4) | (dat1[:,:choicenum3].==3) | (dat1[:,:choicenum3].==4)) & (dat1[:, :level1_hospitals0] .== Lev1 ) & (dat1[:, :level2solo_hospitals0] .== Lev2) & (dat1[:, :level3_hospitals0] .== Lev3))

# will need to duplicate the above separately for entrants.

        # Picks out choices made by firms in that market configuration
        #  choices =  hcat( dat1[all_ind, :choicenum0], dat1[all_ind, :pr_ch_0],  dat1[all_ind, :choicenum1], dat1[all_ind, :pr_ch_1], dat1[all_ind, :choicenum2], dat1[all_ind, :pr_ch_2], dat1[all_ind, :choicenum3], dat1[all_ind, :pr_ch_3], dat1[ all_ind , :act_int], dat1[ all_ind , :act_solo])


        # Identify the Level 1's, level 2's and level 3's
        # lev3 = choices[:, 9].== 1;
        # lev2 = choices[:, 10].== 1;
        # lev1 = ((choices[:,9].==0).*(choices[:,10].==0));
    # Level 1
        if (Lev1 > 0)
          choices1 =  hcat( dat1[ind1, :choicenum0], dat1[ind1, :pr_ch_0],  dat1[ind1, :choicenum1], dat1[ind1, :pr_ch_1], dat1[ind1, :choicenum2], dat1[ind1, :pr_ch_2], dat1[ind1, :choicenum3], dat1[ind1, :pr_ch_3], dat1[ind1, :act_int], dat1[ind1, :act_solo])
          if sum(ind1) > 1
            loc1 = choices1[1, :]
            p11  = loc1[findfirst(loc1, 10)+1]
            p12  = loc1[findfirst(loc1, 2 )+1]
            p13  = loc1[findfirst(loc1, 1 )+1]
            p1ex = loc1[findfirst(loc1, 11)+1]
          elseif sum(ind1) == 1
             p11 = choices1[findfirst(choices1, 10)+1]
             p12 = choices1[findfirst(choices1, 2)+1]
             p13 = choices1[findfirst(choices1, 1)+1]
             p1ex = choices1[findfirst(choices1, 11)+1]
           elseif sum(ind1) == 0
             return ProbException
          end
        elseif  (Lev1 == 0) # the condition I want to catch is when Lev1 > 0 but sum ind == 0
          p11 = 0;
          p12 = 0;
          p13 = 0;
          p1ex = 0;
        end
    # Level 2
        if (Lev2>0)
          choices2 =  hcat( dat1[ind2, :choicenum0], dat1[ind2, :pr_ch_0],  dat1[ind2, :choicenum1], dat1[ind2, :pr_ch_1], dat1[ind2, :choicenum2], dat1[ind2, :pr_ch_2], dat1[ind2, :choicenum3], dat1[ind2, :pr_ch_3], dat1[ind2, :act_int], dat1[ind2, :act_solo])
          if sum(ind2) > 1
            loc2 = choices2[1, :] # Takes the first row - I don't know when this might screw up, but think carefully about it.
            p21 = loc2[findfirst(loc2, 5)+1]
            p22 = loc2[findfirst(loc2, 10)+1]
            p23 = loc2[findfirst(loc2, 6)+1]
            p2ex = loc2[findfirst(loc2, 11)+1]
          elseif sum(ind2) == 1
            p21 = choices2[findfirst(choices2, 5)+1]
            p22 = choices2[findfirst(choices2, 10)+1]
            p23 = choices2[findfirst(choices2, 6)+1]
            p2ex = choices2[findfirst(choices2, 11)+1]
          elseif sum(ind2) == 0
            return ProbException
          end
        elseif Lev2 == 0
          p21 = 0;
          p22 = 0;
          p23 = 0;
          p2ex = 0;
        end
  # Level 3
        if (Lev3 > 0)
          choices3 =  hcat( dat1[ind3, :choicenum0], dat1[ind3, :pr_ch_0],  dat1[ind3, :choicenum1], dat1[ind3, :pr_ch_1], dat1[ind3, :choicenum2], dat1[ind3, :pr_ch_2], dat1[ind3, :choicenum3], dat1[ind3, :pr_ch_3], dat1[ind3, :act_int], dat1[ind3, :act_solo])
          if sum(ind3) > 1
            loc3 = choices3[1,:]
            p31  = loc3[findfirst(loc3, 4)+1]
            p32  = loc3[findfirst(loc3, 3)+1]
            p33  = loc3[findfirst(loc3, 10)+1]
            p3ex = loc3[findfirst(loc3, 11)+1]
          elseif sum(ind3) == 1
            p31 = choices[lev3,:][findfirst(choices[lev3, :], 4)+1]
            p32 = choices[lev3,:][findfirst(choices[lev3, :], 3)+1]
            p33 = choices[lev3,:][findfirst(choices[lev3, :], 10)+1]
            p3ex = choices[lev3,:][findfirst(choices[lev3, :], 11)+1]
          elseif sum(ind3) == 0
            return ProbException
          end
        elseif Lev3 == 0
          p31 = 0;
          p32 = 0;
          p33 = 0;
          p3ex = 0;
        end
      end
    return [p11, p12, p13, p1ex, p21, p22, p23, p2ex, p31, p32, p33, p3ex]
end
