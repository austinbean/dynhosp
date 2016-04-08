# Logit Probs Estimator.
regcoeffs = readtable("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Choice Model.csv", header = true);


CHECK THAT ALL OF THESE ARE NUMBERED CORRECTLY!

# Row Locations - enter by hand?
# Choice specific
NeoIntensive = regcoeffs[1, 2]
SoloIntermediate = regcoeffs[2, 2]
level3_hospitals = regcoeffs[3, 2]
level2solo_hospitals = regcoeffs[4, 2]
level1_hospitals = regcoeffs[5, 2]
plevel3_2 = regcoeffs[6, 2]
plevel3_3 = regcoeffs[7, 2]
plevel3_4 = regcoeffs[8, 2]
plevel3_5 = regcoeffs[9, 2]
plevel3_6 = regcoeffs[10, 2]
plevel3_7 = regcoeffs[11, 2]
plevel3_8 = regcoeffs[12, 2]
plevel3_9 = regcoeffs[13, 2]
plevel3_10 = regcoeffs[14, 2]
plevel2_2 = regcoeffs[15, 2]
plevel2_3 = regcoeffs[16, 2]
plevel2_4 = regcoeffs[17, 2]
plevel2_5 = regcoeffs[18, 2]
plevel2_6 = regcoeffs[19, 2]
plevel2_7 = regcoeffs[20, 2]
plevel2_8 = regcoeffs[21, 2]
plevel2_9 = regcoeffs[22, 2]
plevel1_2 = regcoeffs[23, 2]
plevel1_3 = regcoeffs[24, 2]
plevel1_4 = regcoeffs[25, 2]
plevel1_5 = regcoeffs[26, 2]
plevel1_6 = regcoeffs[27, 2]
plevel1_7 = regcoeffs[28, 2]
plevel1_8 = regcoeffs[29, 2]
plevel1_9 = regcoeffs[30, 2]
plevel1_10 = regcoeffs[31, 2]
# Casevars 13
lev105_13 = regcoeffs[32, 2]
lev205_13 = regcoeffs[33, 2]
lev305_13 = regcoeffs[34, 2]
lev1515_13 = regcoeffs[35, 2]
lev2515_13 = regcoeffs[36, 2]
lev3515_13 = regcoeffs[37, 2]
lev11525_13 = regcoeffs[38, 2]
lev21525_13 = regcoeffs[39, 2]
lev31525_13 = regcoeffs[40, 2]
cons_13 = regcoeffs[41, 2]
# Casevars 21
lev105_21 = regcoeffs[42, 2]
lev205_21 = regcoeffs[43, 2]
lev305_21 = regcoeffs[44, 2]
lev1515_21 = regcoeffs[45, 2]
lev2515_21 = regcoeffs[46, 2]
lev3515_21 = regcoeffs[47, 2]
lev11525_21 = regcoeffs[48, 2]
lev21525_21 = regcoeffs[49, 2]
lev31525_21 = regcoeffs[50, 2]
cons_21 = regcoeffs[51, 2]
# Casevars 23
lev105_23 = regcoeffs[52, 2]
lev205_23 = regcoeffs[53, 2]
lev305_23 = regcoeffs[54, 2]
lev1515_23 = regcoeffs[55, 2]
lev2515_23 = regcoeffs[56, 2]
lev3515_23 = regcoeffs[57, 2]
lev11525_23 = regcoeffs[58, 2]
lev21525_23 = regcoeffs[59, 2]
lev31525_23 = regcoeffs[60, 2]
cons_23 = regcoeffs[61, 2]
# Casevars 13
lev105_13 = regcoeffs[62, 2]
lev205_13 = regcoeffs[63, 2]
lev305_13 = regcoeffs[64, 2]
lev1515_13 = regcoeffs[65, 2]
lev2515_13 = regcoeffs[66, 2]
lev3515_13 = regcoeffs[67, 2]
lev11525_13 = regcoeffs[68, 2]
lev21525_13 = regcoeffs[69, 2]
lev31525_13 = regcoeffs[70, 2]
cons_13 = regcoeffs[71, 2]
# Casevars 12
lev105_12 = regcoeffs[72, 2]
lev205_12 = regcoeffs[73, 2]
lev305_12 = regcoeffs[74, 2]
lev1515_12 = regcoeffs[75, 2]
lev2515_12 = regcoeffs[76, 2]
lev3515_12 = regcoeffs[77, 2]
lev11525_12 = regcoeffs[78, 2]
lev21525_12 = regcoeffs[79, 2]
lev31525_12 = regcoeffs[80, 2]
cons_12 = regcoeffs[81, 2]
# Casevars 32
lev105_32 = regcoeffs[82, 2]
lev205_32 = regcoeffs[83, 2]
lev305_32 = regcoeffs[84, 2]
lev1515_32 = regcoeffs[85, 2]
lev2515_32 = regcoeffs[86, 2]
lev3515_32 = regcoeffs[87, 2]
lev11525_32 = regcoeffs[88, 2]
lev21525_32 = regcoeffs[89, 2]
lev31525_32 = regcoeffs[90, 2]
cons_32 = regcoeffs[91, 2]
# Casevars 31
lev105_31 = regcoeffs[92, 2]
lev205_31 = regcoeffs[93, 2]
lev305_31 = regcoeffs[94, 2]
lev1515_31 = regcoeffs[95, 2]
lev2515_31 = regcoeffs[96, 2]
lev3515_31 = regcoeffs[97, 2]
lev11525_31 = regcoeffs[98, 2]
lev21525_31 = regcoeffs[99, 2]
lev31525_31 = regcoeffs[100, 2]
cons_31 = regcoeffs[101, 2]

# Basic Elements:
NeoIntensive SoloIntermediate level3_hospitals level2solo_hospitals level1_hospitals
# Polynomial Level 3:
plevel3_2 plevel3_3 plevel3_4 plevel3_5 plevel3_6 plevel3_7 plevel3_8 plevel3_9 plevel3_10
# Polynomial Level 2:
plevel2_2 plevel2_3 plevel2_4 plevel2_5 plevel2_6 plevel2_7 plevel2_8 plevel2_9
# Polynomial Level 1:
plevel1_2 plevel1_3 plevel1_4 plevel1_5 plevel1_6 plevel1_7 plevel1_8 plevel1_9 plevel1_10

# Case vars
case13 = [lev105_13 lev205_13 lev305_13 lev1515_13 lev2515_13 lev3515_13 lev11525_13 lev21525_13 lev31525_13 cons_13 ]
case21 = [lev105_21 lev205_21 lev305_21 lev1515_21 lev2515_21 lev3515_21 lev11525_21 lev21525_21 lev31525_21 cons_21 ]
case23 = [lev105_23 lev205_23 lev305_23 lev1515_23 lev2515_23 lev3515_23 lev11525_23 lev21525_23 lev31525_23 cons_23 ]
case13 = [lev105_13 lev205_13 lev305_13 lev1515_13 lev2515_13 lev3515_13 lev11525_13 lev21525_13 lev31525_13 cons_13 ]
[lev105_12 lev205_12 lev305_12 lev1515_12 lev2515_12 lev3515_12 lev11525_12 lev21525_12 lev31525_12 cons_12 ]
[lev105_32 lev205_32 lev305_32 lev1515_32 lev2515_32 lev3515_32 lev11525_32 lev21525_32 lev31525_32 cons_32 ]
[lev105_31 lev205_31 lev305_31 lev1515_31 lev2515_31 lev3515_31 lev11525_31 lev21525_31 lev31525_31 cons_31 ]



#=

# These just write the names of the variables above out in an easy way - no useful functionality

for row in 1:size(regcoeffs)[1]
  startstr = regcoeffs[row,:var][1:search(regcoeffs[row,:var], ':')-1]
  if startstr == "choicenum"
    println( regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], " = regcoeffs[$row, 2]")
  elseif startstr == "To_3_from_1"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_13"]), " = regcoeffs[$row, 2]")
  elseif startstr == "To_2_from_1"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_12" ]), " = regcoeffs[$row, 2]")
  elseif startstr == "To_2_from_3"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_32"]), " = regcoeffs[$row, 2]")
  elseif startstr == "To_1_from_3"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_31"]), " = regcoeffs[$row, 2]")
  elseif startstr == "To_1_from_2"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_21"]), " = regcoeffs[$row, 2]")
  elseif startstr == "To_3_from_2"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_23"]), " = regcoeffs[$row, 2]")
  elseif startstr == "Exit"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_EX"]), " = regcoeffs[$row, 2]")
  end
end


for row in 1:size(regcoeffs)[1]
  startstr = regcoeffs[row,:var][1:search(regcoeffs[row,:var], ':')-1]
  if startstr == "choicenum"
    println( regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])])
  elseif startstr == "To_3_from_1"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_13 "]))
  elseif startstr == "To_2_from_1"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_12 " ]))
  elseif startstr == "To_2_from_3"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_32 "]))
  elseif startstr == "To_1_from_3"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_31 "]))
  elseif startstr == "To_1_from_2"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_21 "]))
  elseif startstr == "To_3_from_2"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_23 "]))
  elseif startstr == "Exit"
    println(join([regcoeffs[row,:var][search(regcoeffs[row,:var], ':')+1:length(regcoeffs[row,:var])], "_EX "]))
  end
end

=#
