# Logit Probs Estimator.



type ValueException <: Exception end

function states1(lev1::Int64, lev2::Int64, lev3::Int64)
  if !((lev1 >= 0) & (lev2 >= 0) & (lev3 >= 0))
    return ValueException
  end
  # format: Intensive, Intermediate, Lev1, Lev2, Lev3
    # st11 =  [0 0 lev1 lev2 lev3]
    # st12 =  [0 1 max(lev1-1, 0) lev2+1 lev3]
    # st13 =  [1 0 max(lev1-1,0) lev2 lev3+1]
    # st1ex = [0 0 max(lev1-1,0) lev2 lev3]
  Result = [0 0 lev1 lev2 lev3; 0 1 max(lev1-1, 0) lev2+1 lev3; 1 0 max(lev1-1,0) lev2 lev3+1; 0 0 max(lev1-1,0) lev2 lev3]
  return Result
end


function states2(lev1::Int64, lev2::Int64, lev3::Int64)
  if !((lev1 >= 0) & (lev2 >= 0) & (lev3 >= 0))
    return ValueException
  end
  # format: Intensive, Intermediate, Lev1, Lev2, Lev3
    # st21  = [0 0 lev1+1 max(lev2-1, 0), lev3]
    # st22  = [0 1 lev1 lev2 lev3]
    # st23  = [1 0 lev1 max(lev2-1, 0) lev3]
    # st2ex = [0 0 lev1 max(lev2-1, 0) lev3]
  Result = [0 0 lev1+1 max(lev2-1, 0) lev3; 0 1 lev1 lev2 lev3; 1 0 lev1 max(lev2-1, 0) lev3; 0 0 lev1 max(lev2-1, 0) lev3]
  return Result
end


function states3(lev1::Int64, lev2::Int64, lev3::Int64)
  if !((lev1 >= 0) & (lev2 >= 0) & (lev3 >= 0))
    return ValueException
  end
  # format: Intensive, Intermediate, Lev1, Lev2, Lev3
    # st31  = [0 0 lev1+1 lev2 max(lev3-1,0)]
    # st32  = [0 1 lev1 lev2+1 max(lev3-1,0)]
    # st33  = [1 0 lev1 lev2 lev3]
    # st3ex = [0 0 lev1 lev2 max(lev3-1,0)]
  Result = [0 0 lev1+1 lev2 max(lev3-1,0); 0 1 lev1 lev2+1 max(lev3-1,0); 1 0 lev1 lev2 lev3; 0 0 lev1 lev2 max(lev3-1,0)]
  return Result
end


function poly( lev, coeffs)
  degree = maximum(size(coeffs))+1
  return coeffs*(lev.^collect(2:degree))
end


function NewPoly(numlev::Int64, level::Int64)  # first arg is number of facs at whatever level, second is level.
  # poly1 = [plevel1_2 plevel1_3 plevel1_4 plevel1_5 plevel1_6 plevel1_7 plevel1_8 plevel1_9 plevel1_10]
  # poly2 = [plevel2_2 plevel2_3 plevel2_4 plevel2_5 plevel2_6 plevel2_7 plevel2_8 plevel2_9]
  # poly3 = [plevel3_2 plevel3_3 plevel3_4 plevel3_5 plevel3_6 plevel3_7 plevel3_8 plevel3_9 plevel3_10]
  coeffs3::Array{Float64,2} = [2.1085  -0.90306  0.233774  -0.0384249  0.00407243  -0.000278157  1.18195e-5  -2.83e-7  2.9e-9] # ProjectModule.Poly3
  coeffs2::Array{Float64,2} = [10.5749  -9.82811  4.88716  -1.40579  0.238834  -0.0235274  0.00123874  -2.69145e-5] #ProjectModule.Poly2
  coeffs1::Array{Float64,2} = [-0.344434  0.0463555  -0.0371076  0.0159125  -0.00313691  0.000333805  -1.99556e-5  6.31e-7  -8.22e-9] #ProjectModule.Poly1
  outp::Float64 = 0.0
  if level == 1
    for i = 1:maximum(size(coeffs1))
      outp+=coeffs1[i]*numlev^(i+1)
    end 
  elseif level == 2
    for i = 1:maximum(size(coeffs2))
      outp+=coeffs2[i]*numlev^(i+1)
    end 
  elseif level == 3
    for i = 1:maximum(size(coeffs3))
      outp+=coeffs3[i]*numlev^(i+1)
    end     
  end 
  return outp
end 

function BasicProd(states::Array{Int64}) 
  # basic = [NeoIntensive SoloIntermediate  level1_hospitals level2solo_hospitals level3_hospitals]
  coeffs::Array{Float64,2}=[5.52019  3.25336  2.86751  -1.74685  -4.2777] # these are the coeffs in "basic" above.
  outp::Float64 = 0.0
  for i = 1:maximum(size(states))
    outp += states[i]*coeffs[i]
  end 
  return outp
end 

function NeighborsProd(cases::Array{Float64,2}, neighbors::Array{Int64})
  outp::Float64 = 0.0
  for i = 1:maximum(size(cases))
    outp += cases[i]*neighbors[i]
  end 
  return outp
end 



function logitest(ownlev::Tuple, lev1::Int64, lev2::Int64, lev3::Int64, neighbors::Array )
  push!(neighbors,1) # adds a 1 at the end of the neighbors to handle the choice-specific constant.
  # here is where to insert the constants... NeighborsProd uses caseXY and neighbors  
  # Casevars 13
  lev105_13::Float64 = ProjectModule.regcoeffs[32, 2]
  lev205_13::Float64 = ProjectModule.regcoeffs[33, 2]
  lev305_13::Float64 = ProjectModule.regcoeffs[34, 2]
  lev1515_13::Float64 = ProjectModule.regcoeffs[35, 2]
  lev2515_13::Float64 = ProjectModule.regcoeffs[36, 2]
  lev3515_13::Float64 = ProjectModule.regcoeffs[37, 2]
  lev11525_13::Float64 = ProjectModule.regcoeffs[38, 2]
  lev21525_13::Float64 = ProjectModule.regcoeffs[39, 2]
  lev31525_13::Float64 = ProjectModule.regcoeffs[40, 2]
  cons_13::Float64 = ProjectModule.regcoeffs[41, 2]
  # Casevars 12
  lev105_12::Float64 = ProjectModule.regcoeffs[42, 2]
  lev205_12::Float64 = ProjectModule.regcoeffs[43, 2]
  lev305_12::Float64 = ProjectModule.regcoeffs[44, 2]
  lev1515_12::Float64 = ProjectModule.regcoeffs[45, 2]
  lev2515_12::Float64 = ProjectModule.regcoeffs[46, 2]
  lev3515_12::Float64 = ProjectModule.regcoeffs[47, 2]
  lev11525_12::Float64 = ProjectModule.regcoeffs[48, 2]
  lev21525_12::Float64 = ProjectModule.regcoeffs[49, 2]
  lev31525_12::Float64 = ProjectModule.regcoeffs[50, 2]
  cons_12::Float64 = ProjectModule.regcoeffs[51, 2]
  # Casevars 32
  lev105_32::Float64 = ProjectModule.regcoeffs[52, 2]
  lev205_32::Float64 = ProjectModule.regcoeffs[53, 2]
  lev305_32::Float64 = ProjectModule.regcoeffs[54, 2]
  lev1515_32::Float64 = ProjectModule.regcoeffs[55, 2]
  lev2515_32::Float64 = ProjectModule.regcoeffs[56, 2]
  lev3515_32::Float64 = ProjectModule.regcoeffs[57, 2]
  lev11525_32::Float64 = ProjectModule.regcoeffs[58, 2]
  lev21525_32::Float64 = ProjectModule.regcoeffs[59, 2]
  lev31525_32::Float64 = ProjectModule.regcoeffs[60, 2]
  cons_32::Float64 = ProjectModule.regcoeffs[61, 2]
  # Casevars 31
  lev105_31::Float64 = ProjectModule.regcoeffs[62, 2]
  lev205_31::Float64 = ProjectModule.regcoeffs[63, 2]
  lev305_31::Float64 = ProjectModule.regcoeffs[64, 2]
  lev1515_31::Float64 = ProjectModule.regcoeffs[65, 2]
  lev2515_31::Float64 = ProjectModule.regcoeffs[66, 2]
  lev3515_31::Float64 = ProjectModule.regcoeffs[67, 2]
  lev11525_31::Float64 = ProjectModule.regcoeffs[68, 2]
  lev21525_31::Float64 = ProjectModule.regcoeffs[69, 2]
  lev31525_31::Float64 = ProjectModule.regcoeffs[70, 2]
  cons_31::Float64 = ProjectModule.regcoeffs[71, 2]
  # Casevars 21
  lev105_21::Float64 = ProjectModule.regcoeffs[72, 2]
  lev205_21::Float64 = ProjectModule.regcoeffs[73, 2]
  lev305_21::Float64 = ProjectModule.regcoeffs[74, 2]
  lev1515_21::Float64 = ProjectModule.regcoeffs[75, 2]
  lev2515_21::Float64 = ProjectModule.regcoeffs[76, 2]
  lev3515_21::Float64 = ProjectModule.regcoeffs[77, 2]
  lev11525_21::Float64 = ProjectModule.regcoeffs[78, 2]
  lev21525_21::Float64 = ProjectModule.regcoeffs[79, 2]
  lev31525_21::Float64 = ProjectModule.regcoeffs[80, 2]
  cons_21::Float64 = ProjectModule.regcoeffs[81, 2]
  # Casevars 23
  lev105_23::Float64 = ProjectModule.regcoeffs[82, 2]
  lev205_23::Float64 = ProjectModule.regcoeffs[83, 2]
  lev305_23::Float64 = ProjectModule.regcoeffs[84, 2]
  lev1515_23::Float64 = ProjectModule.regcoeffs[85, 2]
  lev2515_23::Float64 = ProjectModule.regcoeffs[86, 2]
  lev3515_23::Float64 = ProjectModule.regcoeffs[87, 2]
  lev11525_23::Float64 = ProjectModule.regcoeffs[88, 2]
  lev21525_23::Float64 = ProjectModule.regcoeffs[89, 2]
  lev31525_23::Float64 = ProjectModule.regcoeffs[90, 2]
  cons_23::Float64 = ProjectModule.regcoeffs[91, 2]
  # casevars exit
  lev105_EX::Float64 = ProjectModule.regcoeffs[92, 2]
  lev205_EX::Float64 = ProjectModule.regcoeffs[93, 2]
  lev305_EX::Float64 = ProjectModule.regcoeffs[94, 2]
  lev1515_EX::Float64 = ProjectModule.regcoeffs[95, 2]
  lev2515_EX::Float64 = ProjectModule.regcoeffs[96, 2]
  lev3515_EX::Float64 = ProjectModule.regcoeffs[97, 2]
  lev11525_EX::Float64 = ProjectModule.regcoeffs[98, 2]
  lev21525_EX::Float64 = ProjectModule.regcoeffs[99, 2]
  lev31525_EX::Float64 = ProjectModule.regcoeffs[100, 2]
  cons_EX::Float64 = ProjectModule.regcoeffs[101, 2]
  case13::Array{Float64,2} = [lev105_13 lev205_13 lev305_13 lev1515_13 lev2515_13 lev3515_13 lev11525_13 lev21525_13 lev31525_13 cons_13]
  case12::Array{Float64,2} = [lev105_12 lev205_12 lev305_12 lev1515_12 lev2515_12 lev3515_12 lev11525_12 lev21525_12 lev31525_12 cons_12]
  case32::Array{Float64,2} = [lev105_32 lev205_32 lev305_32 lev1515_32 lev2515_32 lev3515_32 lev11525_32 lev21525_32 lev31525_32 cons_32]
  case31::Array{Float64,2} = [lev105_31 lev205_31 lev305_31 lev1515_31 lev2515_31 lev3515_31 lev11525_31 lev21525_31 lev31525_31 cons_31]
  case21::Array{Float64,2} = [lev105_21 lev205_21 lev305_21 lev1515_21 lev2515_21 lev3515_21 lev11525_21 lev21525_21 lev31525_21 cons_21]
  case23::Array{Float64,2} = [lev105_23 lev205_23 lev305_23 lev1515_23 lev2515_23 lev3515_23 lev11525_23 lev21525_23 lev31525_23 cons_23]
  caseEX::Array{Float64,2} = [lev105_EX lev205_EX lev305_EX lev1515_EX lev2515_EX lev3515_EX lev11525_EX lev21525_EX lev31525_EX cons_EX]
  if ownlev == (0,0) #level 1
    # choices - continue, 12, 13, ex
    retst1 = states1(lev1, lev2, lev3)  
    c11 = BasicProd(retst1[1,:])+NewPoly(retst1[1, 3], 1) + NewPoly(retst1[1, 4], 2) + NewPoly(retst1[1, 5], 3)
    c12 = BasicProd(retst1[2,:])+NewPoly(retst1[2, 3], 1) + NewPoly(retst1[2, 4], 2) + NewPoly(retst1[2, 5], 3) + NeighborsProd(case12,neighbors)
    c13 = BasicProd(retst1[3,:])+NewPoly(retst1[3, 3], 1) + NewPoly(retst1[3, 4], 2) + NewPoly(retst1[3, 5], 3) + NeighborsProd(case13,neighbors)
    c1ex = BasicProd(retst1[4,:])+NewPoly(retst1[4, 3], 1) + NewPoly(retst1[4, 4], 2) + NewPoly(retst1[4, 5], 3) + NeighborsProd(caseEX,neighbors)
    # now compute all of these, but keep in mind that they are differenced.
    # dots don't actually make a difference here.
    p12  = exp.(c12 - c11)./(1.+exp.(c12-c11).+ exp.(c13 - c11).+exp.(c1ex - c11))
    p13  = exp.(c13 - c11)./(1.+exp.(c12-c11).+exp.(c13-c11).+exp.(c1ex - c11))
    p1ex = exp.(c1ex - c11)./(1.+exp.(c12-c11).+exp.(c13-c11).+exp.(c1ex - c11))
    p11 = 1 - p12 - p13 - p1ex
    return [p11 p12 p13 p1ex]
  elseif ownlev == (1,0) # level 2
    # choices - 21, continue, 23, ex
    retst2 = states2(lev1, lev2, lev3) 
    c21 = BasicProd(retst2[1,:])+NewPoly(retst2[1,3],1) + NewPoly(retst2[1,4], 2) + NewPoly(retst2[1, 5], 3) + NeighborsProd(case21,neighbors)
    c22 = BasicProd(retst2[2,:])+NewPoly(retst2[2,3],1) + NewPoly(retst2[2,4], 2) + NewPoly(retst2[2, 5], 3)
    c23 = BasicProd(retst2[3,:])+NewPoly(retst2[3,3],1) + NewPoly(retst2[3,4], 2) + NewPoly(retst2[3, 5], 3) + NeighborsProd(case23,neighbors)
    c2ex = BasicProd(retst2[4,:])+NewPoly(retst2[4,3],1) + NewPoly(retst2[4,4], 2) + NewPoly(retst2[4, 5], 3) + NeighborsProd(caseEX,neighbors)
    # compute all of them:
    p21 = exp(c21 - c22)/(1 + exp(c21 - c22) + exp(c23 - c22) + exp(c2ex - c22))
    p23 = exp(c23 - c22)/(1 + exp(c21 - c22) + exp(c23 - c22) + exp(c2ex - c22))
    p2ex = exp(c2ex - c22)/(1 + exp(c21-c22) + exp(c23 - c22) + exp(c2ex - c22))
    p22 = 1 - p21 - p23 - p2ex
    return [p21 p22 p23 p2ex]
  elseif ownlev == (0,1) # level 3
    # choices - 31, 32, continue, ex
    retst3 = states3(lev1, lev2, lev3) 
    c31 = BasicProd(retst3[1,:])+NewPoly(retst3[1, 3], 1) + NewPoly(retst3[1, 4], 2) + NewPoly(retst3[1, 5], 3) + NeighborsProd(case31,neighbors)
    c32 = BasicProd(retst3[2,:])+NewPoly(retst3[2, 3], 1) + NewPoly(retst3[2, 4], 2) + NewPoly(retst3[2, 5], 3) + NeighborsProd(case32,neighbors)
    c33 = BasicProd(retst3[3,:])+NewPoly(retst3[3, 3], 1) + NewPoly(retst3[3, 4], 2) + NewPoly(retst3[3, 5], 3)
    c3ex = BasicProd(retst3[4,:])+ NewPoly(retst3[4, 3], 1) + NewPoly(retst3[4, 4], 2) + NewPoly(retst3[4, 5], 3) + NeighborsProd(caseEX,neighbors)
    # compute:
    p31 = exp(c31 - c33)/(1 + exp(c31 - c33) + exp(c32 - c33) + exp(c3ex - c33))
    p32 = exp(c32 - c33)/(1 + exp(c31 - c33) + exp(c32 - c33) + exp(c3ex - c33))
    p3ex = exp(c3ex - c33)/(1 + exp(c31 - c33) + exp(c32 - c33) + exp(c3ex - c33))
    p33 = 1 - p31 - p32 - p3ex
    return [p31 p32 p33 p3ex]
  else
    return ValueException
  end
end



#=
# Row Locations - enter by hand?
# Choice specific
NeoIntensive = ProjectModule.regcoeffs[1, 2]
SoloIntermediate = ProjectModule.regcoeffs[2, 2]
level3_hospitals = ProjectModule.regcoeffs[3, 2]
level2solo_hospitals = ProjectModule.regcoeffs[4, 2]
level1_hospitals = ProjectModule.regcoeffs[5, 2]
=#

#=
plevel3_2 = ProjectModule.regcoeffs[6, 2]
plevel3_3 = ProjectModule.regcoeffs[7, 2]
plevel3_4 = ProjectModule.regcoeffs[8, 2]
plevel3_5 = ProjectModule.regcoeffs[9, 2]
plevel3_6 = ProjectModule.regcoeffs[10, 2]
plevel3_7 = ProjectModule.regcoeffs[11, 2]
plevel3_8 = ProjectModule.regcoeffs[12, 2]
plevel3_9 = ProjectModule.regcoeffs[13, 2]
plevel3_10 = ProjectModule.regcoeffs[14, 2]
plevel2_2 = ProjectModule.regcoeffs[15, 2]
plevel2_3 = ProjectModule.regcoeffs[16, 2]
plevel2_4 = ProjectModule.regcoeffs[17, 2]
plevel2_5 = ProjectModule.regcoeffs[18, 2]
plevel2_6 = ProjectModule.regcoeffs[19, 2]
plevel2_7 = ProjectModule.regcoeffs[20, 2]
plevel2_8 = ProjectModule.regcoeffs[21, 2]
plevel2_9 = ProjectModule.regcoeffs[22, 2]
plevel1_2 = ProjectModule.regcoeffs[23, 2]
plevel1_3 = ProjectModule.regcoeffs[24, 2]
plevel1_4 = ProjectModule.regcoeffs[25, 2]
plevel1_5 = ProjectModule.regcoeffs[26, 2]
plevel1_6 = ProjectModule.regcoeffs[27, 2]
plevel1_7 = ProjectModule.regcoeffs[28, 2]
plevel1_8 = ProjectModule.regcoeffs[29, 2]
plevel1_9 = ProjectModule.regcoeffs[30, 2]
plevel1_10 = ProjectModule.regcoeffs[31, 2]
=#

#=
# Casevars 13
lev105_13 = ProjectModule.regcoeffs[32, 2]
lev205_13 = ProjectModule.regcoeffs[33, 2]
lev305_13 = ProjectModule.regcoeffs[34, 2]
lev1515_13 = ProjectModule.regcoeffs[35, 2]
lev2515_13 = ProjectModule.regcoeffs[36, 2]
lev3515_13 = ProjectModule.regcoeffs[37, 2]
lev11525_13 = ProjectModule.regcoeffs[38, 2]
lev21525_13 = ProjectModule.regcoeffs[39, 2]
lev31525_13 = ProjectModule.regcoeffs[40, 2]
cons_13 = ProjectModule.regcoeffs[41, 2]
# Casevars 12
lev105_12 = ProjectModule.regcoeffs[42, 2]
lev205_12 = ProjectModule.regcoeffs[43, 2]
lev305_12 = ProjectModule.regcoeffs[44, 2]
lev1515_12 = ProjectModule.regcoeffs[45, 2]
lev2515_12 = ProjectModule.regcoeffs[46, 2]
lev3515_12 = ProjectModule.regcoeffs[47, 2]
lev11525_12 = ProjectModule.regcoeffs[48, 2]
lev21525_12 = ProjectModule.regcoeffs[49, 2]
lev31525_12 = ProjectModule.regcoeffs[50, 2]
cons_12 = ProjectModule.regcoeffs[51, 2]
# Casevars 32
lev105_32 = ProjectModule.regcoeffs[52, 2]
lev205_32 = ProjectModule.regcoeffs[53, 2]
lev305_32 = ProjectModule.regcoeffs[54, 2]
lev1515_32 = ProjectModule.regcoeffs[55, 2]
lev2515_32 = ProjectModule.regcoeffs[56, 2]
lev3515_32 = ProjectModule.regcoeffs[57, 2]
lev11525_32 = ProjectModule.regcoeffs[58, 2]
lev21525_32 = ProjectModule.regcoeffs[59, 2]
lev31525_32 = ProjectModule.regcoeffs[60, 2]
cons_32 = ProjectModule.regcoeffs[61, 2]
# Casevars 31
lev105_31 = ProjectModule.regcoeffs[62, 2]
lev205_31 = ProjectModule.regcoeffs[63, 2]
lev305_31 = ProjectModule.regcoeffs[64, 2]
lev1515_31 = ProjectModule.regcoeffs[65, 2]
lev2515_31 = ProjectModule.regcoeffs[66, 2]
lev3515_31 = ProjectModule.regcoeffs[67, 2]
lev11525_31 = ProjectModule.regcoeffs[68, 2]
lev21525_31 = ProjectModule.regcoeffs[69, 2]
lev31525_31 = ProjectModule.regcoeffs[70, 2]
cons_31 = ProjectModule.regcoeffs[71, 2]
# Casevars 21
lev105_21 = ProjectModule.regcoeffs[72, 2]
lev205_21 = ProjectModule.regcoeffs[73, 2]
lev305_21 = ProjectModule.regcoeffs[74, 2]
lev1515_21 = ProjectModule.regcoeffs[75, 2]
lev2515_21 = ProjectModule.regcoeffs[76, 2]
lev3515_21 = ProjectModule.regcoeffs[77, 2]
lev11525_21 = ProjectModule.regcoeffs[78, 2]
lev21525_21 = ProjectModule.regcoeffs[79, 2]
lev31525_21 = ProjectModule.regcoeffs[80, 2]
cons_21 = ProjectModule.regcoeffs[81, 2]
# Casevars 23
lev105_23 = ProjectModule.regcoeffs[82, 2]
lev205_23 = ProjectModule.regcoeffs[83, 2]
lev305_23 = ProjectModule.regcoeffs[84, 2]
lev1515_23 = ProjectModule.regcoeffs[85, 2]
lev2515_23 = ProjectModule.regcoeffs[86, 2]
lev3515_23 = ProjectModule.regcoeffs[87, 2]
lev11525_23 = ProjectModule.regcoeffs[88, 2]
lev21525_23 = ProjectModule.regcoeffs[89, 2]
lev31525_23 = ProjectModule.regcoeffs[90, 2]
cons_23 = ProjectModule.regcoeffs[91, 2]
# casevars exit
lev105_EX = ProjectModule.regcoeffs[92, 2]
lev205_EX = ProjectModule.regcoeffs[93, 2]
lev305_EX = ProjectModule.regcoeffs[94, 2]
lev1515_EX = ProjectModule.regcoeffs[95, 2]
lev2515_EX = ProjectModule.regcoeffs[96, 2]
lev3515_EX = ProjectModule.regcoeffs[97, 2]
lev11525_EX = ProjectModule.regcoeffs[98, 2]
lev21525_EX = ProjectModule.regcoeffs[99, 2]
lev31525_EX = ProjectModule.regcoeffs[100, 2]
cons_EX = ProjectModule.regcoeffs[101, 2]
=#


#=
# Basic Elements:
basic = [NeoIntensive SoloIntermediate  level1_hospitals level2solo_hospitals level3_hospitals]
# Polynomial Level 3:
poly3 = [plevel3_2 plevel3_3 plevel3_4 plevel3_5 plevel3_6 plevel3_7 plevel3_8 plevel3_9 plevel3_10]
# Polynomial Level 2:
poly2 = [plevel2_2 plevel2_3 plevel2_4 plevel2_5 plevel2_6 plevel2_7 plevel2_8 plevel2_9]
# Polynomial Level 1:
poly1 = [plevel1_2 plevel1_3 plevel1_4 plevel1_5 plevel1_6 plevel1_7 plevel1_8 plevel1_9 plevel1_10]
=#


# Case vars
#=
case13 = [lev105_13 lev205_13 lev305_13 lev1515_13 lev2515_13 lev3515_13 lev11525_13 lev21525_13 lev31525_13 cons_13]
case12 = [lev105_12 lev205_12 lev305_12 lev1515_12 lev2515_12 lev3515_12 lev11525_12 lev21525_12 lev31525_12 cons_12]
case32 = [lev105_32 lev205_32 lev305_32 lev1515_32 lev2515_32 lev3515_32 lev11525_32 lev21525_32 lev31525_32 cons_32]
case31 = [lev105_31 lev205_31 lev305_31 lev1515_31 lev2515_31 lev3515_31 lev11525_31 lev21525_31 lev31525_31 cons_31]
case21 = [lev105_21 lev205_21 lev305_21 lev1515_21 lev2515_21 lev3515_21 lev11525_21 lev21525_21 lev31525_21 cons_21]
case23 = [lev105_23 lev205_23 lev305_23 lev1515_23 lev2515_23 lev3515_23 lev11525_23 lev21525_23 lev31525_23 cons_23]
caseEX = [lev105_EX lev205_EX lev305_EX lev1515_EX lev2515_EX lev3515_EX lev11525_EX lev21525_EX lev31525_EX cons_EX]
=#




#=
# Testing the function:
  # Test 1:
  logitest(( dataf[1,:act_int], dataf[1,:act_solo]), dataf[1,:level1_hospitals0], dataf[1,:level2solo_hospitals0], dataf[1,:level3_hospitals0],[dataf[1,:lev105], dataf[1,:lev205], dataf[1,:lev305], dataf[1,:lev1515], dataf[1,:lev2515], dataf[1,:lev3515], dataf[1,:lev11525], dataf[1,:lev21525], dataf[1,:lev31525]] )
  # check against:
  [dataf[1,:choicenum0], dataf[1,:pr_ch_0], dataf[1,:choicenum1], dataf[1,:pr_ch_1], dataf[1,:choicenum2], dataf[1, :pr_ch_2], dataf[1,:choicenum3], dataf[1,:pr_ch_3]]
  # Test 2:
  logitest(( dataf[170,:act_int], dataf[170,:act_solo]), dataf[170,:level1_hospitals0], dataf[170,:level2solo_hospitals0], dataf[170,:level3_hospitals0],[dataf[170,:lev105], dataf[170,:lev205], dataf[170,:lev305], dataf[170,:lev1515], dataf[170,:lev2515], dataf[170,:lev3515], dataf[170,:lev11525], dataf[170,:lev21525], dataf[170,:lev31525]] )
  # check against:
  [dataf[170,:choicenum0], dataf[170,:pr_ch_0], dataf[170,:choicenum1], dataf[170,:pr_ch_1], dataf[170,:choicenum2], dataf[170, :pr_ch_2], dataf[170,:choicenum3], dataf[170,:pr_ch_3]]
  # Test 3:
  logitest(( dataf[5070,:act_int], dataf[5070,:act_solo]), dataf[5070,:level1_hospitals0], dataf[5070,:level2solo_hospitals0], dataf[5070,:level3_hospitals0],[dataf[5070,:lev105], dataf[5070,:lev205], dataf[5070,:lev305], dataf[5070,:lev1515], dataf[5070,:lev2515], dataf[5070,:lev3515], dataf[5070,:lev11525], dataf[5070,:lev21525], dataf[5070,:lev31525]] )
  # check against:
  [dataf[5070,:choicenum0], dataf[5070,:pr_ch_0], dataf[5070,:choicenum1], dataf[5070,:pr_ch_1], dataf[5070,:choicenum2], dataf[5070, :pr_ch_2], dataf[5070,:choicenum3], dataf[5070,:pr_ch_3]]
=#

#=
# These just write the names of the variables above out in an easy way - no useful functionality

for row in 1:size(ProjectModule.regcoeffs)[1]
  startstr = ProjectModule.regcoeffs[row,:var][1:search(ProjectModule.regcoeffs[row,:var], ':')-1]
  if startstr == "choicenum"
    println( ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], " = ProjectModule.regcoeffs[$row, 2]")
  elseif startstr == "To_3_from_1"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_13"]), " = ProjectModule.regcoeffs[$row, 2]")
  elseif startstr == "To_2_from_1"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_12" ]), " = ProjectModule.regcoeffs[$row, 2]")
  elseif startstr == "To_2_from_3"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_32"]), " = ProjectModule.regcoeffs[$row, 2]")
  elseif startstr == "To_1_from_3"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_31"]), " = ProjectModule.regcoeffs[$row, 2]")
  elseif startstr == "To_1_from_2"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_21"]), " = ProjectModule.regcoeffs[$row, 2]")
  elseif startstr == "To_3_from_2"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_23"]), " = ProjectModule.regcoeffs[$row, 2]")
  elseif startstr == "Exit"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_EX"]), " = ProjectModule.regcoeffs[$row, 2]")
  end
end


for row in 1:size(ProjectModule.regcoeffs)[1]
  startstr = ProjectModule.regcoeffs[row,:var][1:search(ProjectModule.regcoeffs[row,:var], ':')-1]
  if startstr == "choicenum"
    println( ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])])
  elseif startstr == "To_3_from_1"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_13 "]))
  elseif startstr == "To_2_from_1"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_12 " ]))
  elseif startstr == "To_2_from_3"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_32 "]))
  elseif startstr == "To_1_from_3"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_31 "]))
  elseif startstr == "To_1_from_2"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_21 "]))
  elseif startstr == "To_3_from_2"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_23 "]))
  elseif startstr == "Exit"
    println(join([ProjectModule.regcoeffs[row,:var][search(ProjectModule.regcoeffs[row,:var], ':')+1:length(ProjectModule.regcoeffs[row,:var])], "_EX "]))
  end
end

=#
