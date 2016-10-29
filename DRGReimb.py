

import csv
import numpy as np

# TODO - the reported values are averages: is this getting that right?  Double check!


drgdata = []
with open('/Users/austinbean/Desktop/justreimb.csv', encoding = "ISO-8859-1") as f:
    a = csv.reader(f, delimiter=',')
    for row in a:
        drgdata.append(row)


nams = set()
drgs = set()
sts = set()
for row in drgdata:
    nams.add(row[2])
    drgs.add(row[9])
    sts.add(row[6])

#drgs = 385, 386, 387, 388, 389, 390, 391
#drgs = 789, 790, 791, 792, 793, 794, 795
# TODO - add the DRG's for the mothers as well.  Figure out how much extra the hospital is getting there.

pats2008 =[0, 0, 0, 0, 0, 0, 0]
chrg2008 =[0, 0, 0, 0, 0, 0, 0]
pats2009 = [0, 0, 0, 0, 0, 0, 0]
chrg2009 = [0, 0, 0, 0, 0, 0, 0]
pats2010 = [0, 0, 0, 0, 0, 0, 0]
chrg2010 = [0, 0, 0, 0, 0, 0, 0]
pats2011 = [0, 0, 0, 0, 0, 0, 0]
chrg2011 = [0, 0, 0, 0, 0, 0, 0]
pats2012 = [0, 0, 0, 0, 0, 0, 0]
chrg2012 = [0, 0, 0, 0, 0, 0, 0]


# Create vectors to hold the reimbursement per-person for all year-drg combinations.  Compute mean and SD later.
yr2008d789 = []
yr2009d789 = []
yr2010d789 = []
yr2011d789 = []
yr2012d789 = []

yr2008d790 = []
yr2009d790 = []
yr2010d790 = []
yr2011d790 = []
yr2012d790 = []

yr2008d791 = []
yr2009d791 = []
yr2010d791 = []
yr2011d791 = []
yr2012d791 = []

yr2008d792 = []
yr2009d792 = []
yr2010d792 = []
yr2011d792 = []
yr2012d792 = []

yr2008d793 = []
yr2009d793 = []
yr2010d793 = []
yr2011d793 = []
yr2012d793 = []

yr2008d794 = []
yr2009d794 = []
yr2010d794 = []
yr2011d794 = []
yr2012d794 = []

yr2008d795 = []
yr2009d795 = []
yr2010d795 = []
yr2011d795 = []
yr2012d795 = []



# How many patients are there?
# Just check 2012 for now since I'll report that
# This makes different assumptions for the censored and uncensored values.  cell sizes below 5 are censored, so they can be either 1 ("lo") or 5 ("hi")

yr2012drg789hi = 0
yr2012drg789lo = 0
yr2012drg790hi = 0
yr2012drg790lo = 0
yr2012drg791hi = 0
yr2012drg791lo = 0
yr2012drg792hi = 0
yr2012drg792lo = 0
yr2012drg793hi = 0
yr2012drg793lo = 0
yr2012drg794hi = 0
yr2012drg794lo = 0
yr2012drg795hi = 0
yr2012drg795lo = 0

myr2012drg765hi = 0
myr2012drg765lo = 0
myr2012drg766hi = 0
myr2012drg766lo = 0
myr2012drg774hi = 0
myr2012drg774lo = 0
myr2012drg775hi = 0
myr2012drg775lo = 0
myr2012drg767hi = 0
myr2012drg767lo = 0
myr2012drg768hi = 0
myr2012drg768lo = 0



for row in drgdata:
    if (("HOSPITAL" in row[2])|("HEALTH CLINIC" in row[2]))&(row[6] == "TX"):
        if '2012' in row[8]: # only computing a patient count for 2012 for the time being
            if '< 5' in row[14]: # every mention of "< 5" is recorded as "< <space> 5"
                if eval(row[9]) == 789:
                    yr2012drg789hi += 5
                    yr2012drg789lo += 1
                elif eval(row[9]) == 790:
                    yr2012drg790hi += 5
                    yr2012drg790lo += 1
                elif eval(row[9]) == 791:
                    yr2012drg791hi += 5
                    yr2012drg791lo += 1
                elif eval(row[9]) == 792:
                    yr2012drg792hi += 5
                    yr2012drg792lo += 1
                elif eval(row[9]) == 793:
                    yr2012drg793hi += 5
                    yr2012drg793lo += 1
                elif eval(row[9]) == 794:
                    yr2012drg794hi += 5
                    yr2012drg794lo += 1
                elif eval(row[9]) == 795:
                    yr2012drg795hi += 5
                    yr2012drg795lo += 1
                elif eval(row[9]) == 765:
                    myr2012drg765hi += 5
                    myr2012drg765lo += 1
                elif eval(row[9]) == 766:
                    myr2012drg766hi += 5
                    myr2012drg766lo += 1
                elif eval(row[9]) == 774:
                    myr2012drg774hi += 5
                    myr2012drg774lo += 1
                elif eval(row[9]) == 775:
                    myr2012drg775hi += 5
                    myr2012drg775lo += 1
                elif eval(row[9]) == 767:
                    myr2012drg767hi += 5
                    myr2012drg767lo += 1
                elif eval(row[9]) == 768:
                    myr2012drg768hi += 5
                    myr2012drg768lo += 1
            else: # every mention of "< 5" is recorded as "< <space> 5"
                if eval(row[9]) == 789:
                    yr2012drg789hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    yr2012drg789lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 790:
                    yr2012drg790hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    yr2012drg790lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 791:
                    yr2012drg791hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    yr2012drg791lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 792:
                    yr2012drg792hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    yr2012drg792lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 793:
                    yr2012drg793hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    yr2012drg793lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 794:
                    yr2012drg794hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    yr2012drg794lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 795:
                    yr2012drg795hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    yr2012drg795lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 765:
                    myr2012drg765hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    myr2012drg765lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 766:
                    myr2012drg766hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    myr2012drg766lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 774:
                    myr2012drg774hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    myr2012drg774lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 767:
                    myr2012drg767hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    myr2012drg767lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                elif eval(row[9]) == 768:
                    myr2012drg768hi += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())
                    myr2012drg768lo += eval(row[14].replace('<', '').replace('>','').replace(',','').strip())

variants = set()
for row in drgdata:
    if (("HOSPITAL" in row[2])|("HEALTH CLINIC" in row[2]))&(row[6] == "TX"):
        if ('<' in row[14])&('5' in row[14]):
            variants.add(row[14])


for row in drgdata:
    if (("HOSPITAL" in row[2])|("HEALTH CLINIC" in row[2]))&(row[6] == "TX"):
        if '>' in row[14]:
            row[14] = row[14].replace('>', '').strip()
        elif '<' in row[14]:
            if ('5' in row[14]):
                row[14] = '1' # less than 5 censored, so replace with 1.   Could also replace with 5.
            else:
                row[14] = row[14].replace('<', '').strip()
        if eval(row[9]) == 789:
            if '2008' in row[8]:
                pats2008[0] += eval(row[14].replace(',','').strip()) # number of patients
                chrg2008[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip()) # AVERAGE charge for that number
                yr2008d789.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[0] += eval(row[14].replace(',','').strip())
                chrg2009[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d789.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[0] += eval(row[14].replace(',','').strip())
                chrg2010[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d789.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[0] += eval(row[14].replace(',','').strip())
                chrg2011[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d789.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[0] += eval(row[14].replace(',','').strip())
                chrg2012[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d789.append(eval(row[15].replace('$','').replace(',','').strip()))
        elif eval(row[9]) == 790:
            if '2008' in row[8]:
                pats2008[1] += eval(row[14].replace(',','').strip())
                chrg2008[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d790.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[1] += eval(row[14].replace(',','').strip())
                chrg2009[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d790.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[1] += eval(row[14].replace(',','').strip())
                chrg2010[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d790.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[1] += eval(row[14].replace(',','').strip())
                chrg2011[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d790.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[1] += eval(row[14].replace(',','').strip())
                chrg2012[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d790.append(eval(row[15].replace('$','').replace(',','').strip()))
        elif eval(row[9]) == 791:
            if '2008' in row[8]:
                pats2008[2] += eval(row[14].replace(',','').strip())
                chrg2008[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d791.append(eval(row[15].replace('$','').replace(',','')))
            elif '2009' in row[8]:
                pats2009[2] += eval(row[14].replace(',','').strip())
                chrg2009[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d791.append(eval(row[15].replace('$','').replace(',','')))
            elif '2010' in row[8]:
                pats2010[2] += eval(row[14].replace(',','').strip())
                chrg2010[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d791.append(eval(row[15].replace('$','').replace(',','')))
            elif '2011' in row[8]:
                pats2011[2] += eval(row[14].replace(',','').strip())
                chrg2011[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d791.append(eval(row[15].replace('$','').replace(',','')))
            elif '2012' in row[8]:
                pats2012[2] += eval(row[14].replace(',','').strip())
                chrg2012[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d791.append(eval(row[15].replace('$','').replace(',','')))
        elif eval(row[9]) == 792:
            if '2008' in row[8]:
                pats2008[3] += eval(row[14].replace(',','').strip())
                chrg2008[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d792.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[3] += eval(row[14].replace(',','').strip())
                chrg2009[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d792.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[3] += eval(row[14].replace(',','').strip())
                chrg2010[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d792.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[3] += eval(row[14].replace(',','').strip())
                chrg2011[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d792.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[3] += eval(row[14].replace(',','').strip())
                chrg2012[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d792.append(eval(row[15].replace('$','').replace(',','').strip()))
        elif eval(row[9]) == 793:
            if '2008' in row[8]:
                pats2008[4] += eval(row[14].replace(',','').strip())
                chrg2008[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d793.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[4] += eval(row[14].replace(',','').strip())
                chrg2009[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d793.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[4] += eval(row[14].replace(',','').strip())
                chrg2010[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d793.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[4] += eval(row[14].replace(',','').strip())
                chrg2011[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d793.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[4] += eval(row[14].replace(',','').strip())
                chrg2012[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d793.append(eval(row[15].replace('$','').replace(',','').strip()))
        elif eval(row[9]) == 794:
            if '2008' in row[8]:
                pats2008[5] += eval(row[14].replace(',','').strip())
                chrg2008[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d794.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[5] += eval(row[14].replace(',','').strip())
                chrg2009[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d794.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[5] += eval(row[14].replace(',','').strip())
                chrg2010[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d794.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[5] += eval(row[14].replace(',','').strip())
                chrg2011[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d794.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[5] += eval(row[14].replace(',','').strip())
                chrg2012[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d794.append(eval(row[15].replace('$','').replace(',','').strip()))
        elif eval(row[9]) == 795:
            if '2008' in row[8]:
                pats2008[6] += eval(row[14].replace(',','').strip())
                chrg2008[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d795.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[6] += eval(row[14].replace(',','').strip())
                chrg2009[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d795.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[6] += eval(row[14].replace(',','').strip())
                chrg2010[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d795.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[6] += eval(row[14].replace(',','').strip())
                chrg2011[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d795.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[6] += eval(row[14].replace(',','').strip())
                chrg2012[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d795.append(eval(row[15].replace('$','').replace(',','').strip()))



moms2008 =[0, 0, 0, 0, 0, 0]
mchrg2008 =[0, 0, 0, 0, 0, 0]
moms2009 = [0, 0, 0, 0, 0, 0]
mchrg2009 = [0, 0, 0, 0, 0, 0]
moms2010 = [0, 0, 0, 0, 0, 0]
mchrg2010 = [0, 0, 0, 0, 0, 0]
moms2011 = [0, 0, 0, 0, 0, 0]
mchrg2011 = [0, 0, 0, 0, 0, 0]
moms2012 = [0, 0, 0, 0, 0, 0]
mchrg2012 = [0, 0, 0, 0, 0, 0]

myr2008d765 = []
myr2009d765 = []
myr2010d765 = []
myr2011d765 = []
myr2012d765 = []

myr2008d766 = []
myr2009d766 = []
myr2010d766 = []
myr2011d766 = []
myr2012d766 = []

myr2008d767 = []
myr2009d767 = []
myr2010d767 = []
myr2011d767 = []
myr2012d767 = []

myr2008d768 = []
myr2009d768 = []
myr2010d768 = []
myr2011d768 = []
myr2012d768 = []

myr2008d774 = []
myr2009d774 = []
myr2010d774 = []
myr2011d774 = []
myr2012d774 = []

myr2008d775 = []
myr2009d775 = []
myr2010d775 = []
myr2011d775 = []
myr2012d775 = []

for row in drgdata:
    if (("HOSPITAL" in row[2])|("HEALTH CLINIC" in row[2]))&(row[6] == "TX"):
        if '>' in row[14]:
            row[14] = row[14].replace('>', '').strip()
        elif '<' in row[14]:
            if ('5' in row[14]):
                row[14] = '1' # less than 5 censored, so replace with 1.   Could also replace with 5.
            else:
                row[14] = row[14].replace('<', '').strip()
        if eval(row[9]) == 765:
            if '2008' in row[8]:
                moms2008[0] += eval(row[14].replace(',','').strip())
                mchrg2008[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2008d765.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2009' in row[8]:
                moms2009[0] += eval(row[14].replace(',','').strip())
                mchrg2009[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2009d765.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2010' in row[8]:
                moms2010[0] += eval(row[14].replace(',','').strip())
                mchrg2010[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2010d765.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2011' in row[8]:
                moms2011[0] += eval(row[14].replace(',','').strip())
                mchrg2011[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2011d765.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2012' in row[8]:
                moms2012[0] += eval(row[14].replace(',','').strip())
                mchrg2012[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2012d765.append(eval(row[15].replace('$', '').replace(',','').strip()))
        elif eval(row[9]) == 766:
            if '2008' in row[8]:
                moms2008[1] += eval(row[14].replace(',','').strip())
                mchrg2008[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2008d766.append(eval(row[15].replace('$', '').replace(',', '').strip()))
            elif '2009' in row[8]:
                moms2009[1] += eval(row[14].replace(',','').strip())
                mchrg2009[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2009d766.append(eval(row[15].replace('$', '').replace(',', '').strip()))
            elif '2010' in row[8]:
                moms2010[1] += eval(row[14].replace(',','').strip())
                mchrg2010[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2010d766.append(eval(row[15].replace('$', '').replace(',', '').strip()))
            elif '2011' in row[8]:
                moms2011[1] += eval(row[14].replace(',','').strip())
                mchrg2011[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2011d766.append(eval(row[15].replace('$', '').replace(',', '').strip()))
            elif '2012' in row[8]:
                moms2012[1] += eval(row[14].replace(',','').strip())
                mchrg2012[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2012d766.append(eval(row[15].replace('$', '').replace(',', '').strip()))
        elif eval(row[9]) == 767:
            if '2008' in row[8]:
                moms2008[2] += eval(row[14].replace(',','').strip())
                mchrg2008[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2008d767.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2009' in row[8]:
                moms2009[2] += eval(row[14].replace(',','').strip())
                mchrg2009[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2009d767.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2010' in row[8]:
                moms2010[2] += eval(row[14].replace(',','').strip())
                mchrg2010[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2010d767.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2011' in row[8]:
                moms2011[2] += eval(row[14].replace(',','').strip())
                mchrg2011[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2011d767.append(eval(row[15].replace('$', '').replace(',','').strip()))
            elif '2012' in row[8]:
                moms2012[2] += eval(row[14].replace(',','').strip())
                mchrg2012[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2012d767.append(eval(row[15].replace('$', '').replace(',','').strip()))
        elif eval(row[9]) == 768:
            if '2008' in row[8]:
                moms2008[3] += eval(row[14].replace(',','').strip())
                mchrg2008[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2008d768.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                moms2009[3] += eval(row[14].replace(',','').strip())
                mchrg2009[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2009d768.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                moms2010[3] += eval(row[14].replace(',','').strip())
                mchrg2010[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2010d768.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                moms2011[3] += eval(row[14].replace(',','').strip())
                mchrg2011[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2011d768.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                moms2012[3] += eval(row[14].replace(',','').strip())
                mchrg2012[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2012d768.append(eval(row[15].replace('$','').replace(',','').strip()))
        elif eval(row[9]) == 774:
            if '2008' in row[8]:
                moms2008[4] += eval(row[14].replace(',','').strip())
                mchrg2008[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2008d774.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2009' in row[8]:
                moms2009[4] += eval(row[14].replace(',','').strip())
                mchrg2009[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2009d774.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2010' in row[8]:
                moms2010[4] += eval(row[14].replace(',','').strip())
                mchrg2010[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2010d774.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2011' in row[8]:
                moms2011[4] += eval(row[14].replace(',','').strip())
                mchrg2011[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2011d774.append(eval(row[15].replace('$','').replace(',','').strip()))
            elif '2012' in row[8]:
                moms2012[4] += eval(row[14].replace(',','').strip())
                mchrg2012[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2012d774.append(eval(row[15].replace('$','').replace(',','').strip()))
        elif eval(row[9]) == 775:
            if '2008' in row[8]:
                moms2008[5] += eval(row[14].replace(',','').strip())
                mchrg2008[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2008d775.append(eval(row[15].replace('$','').replace(',', '').strip()))
            elif '2009' in row[8]:
                moms2009[5] += eval(row[14].replace(',','').strip())
                mchrg2009[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2009d775.append(eval(row[15].replace('$','').replace(',', '').strip()))
            elif '2010' in row[8]:
                moms2010[5] += eval(row[14].replace(',','').strip())
                mchrg2010[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2010d775.append(eval(row[15].replace('$','').replace(',', '').strip()))
            elif '2011' in row[8]:
                moms2011[5] += eval(row[14].replace(',','').strip())
                mchrg2011[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2011d775.append(eval(row[15].replace('$','').replace(',', '').strip()))
            elif '2012' in row[8]:
                moms2012[5] += eval(row[14].replace(',','').strip())
                mchrg2012[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                myr2012d775.append(eval(row[15].replace('$','').replace(',', '').strip()))


# Prints the results of the upper loop:

labelsc= ["385/789 Neonates died or transferred", "386/790 Extreme immaturity/respiratory distress", "387/791 Prematurity w/ Major Problems", "388/792 Prematurity w/out major problems", "389/793 Full term neonate w/ major problems", "390/794 Neonate w/ other significant problems", "391/795 Normal Newborn"]

avg2008 = [0,0,0,0,0,0,0]
avg2009 = [0,0,0,0,0,0,0]
avg2010 = [0,0,0,0,0,0,0]
avg2011 = [0,0,0,0,0,0,0]
avg2012 = [0,0,0,0,0,0,0]

for i in range(0,7) :
    avg2008[i] = chrg2008[i]/pats2008[i]
    avg2009[i] = chrg2009[i]/pats2009[i]
    avg2010[i] = chrg2010[i]/pats2010[i]
    avg2011[i] = chrg2011[i]/pats2011[i]
    avg2012[i] = chrg2012[i]/pats2012[i]


for i in range(0,7):
    print(2008, labelsc[i], avg2008[i])
    print(2009, labelsc[i], avg2009[i])
    print(2010, labelsc[i], avg2010[i])
    print(2011, labelsc[i], avg2011[i])
    print(2012, labelsc[i], avg2012[i])


print(2008, labelsc[0], "mean", np.mean(yr2008d789), "std", np.std(yr2008d789))
print(2009, labelsc[0], "mean", np.mean(yr2009d789), "std", np.std(yr2009d789))
print(2010, labelsc[0], "mean", np.mean(yr2010d789), "std", np.std(yr2010d789))
print(2011, labelsc[0], "mean", np.mean(yr2011d789), "std", np.std(yr2011d789))
print(2012, labelsc[0], " mean ", np.mean(yr2012d789), " std ", np.std(yr2012d789), " num patients low ", yr2012drg789lo, " num patients high ", yr2012drg789hi)

print(2008, labelsc[1], "mean", np.mean(yr2008d790), "std", np.std(yr2008d790))
print(2009, labelsc[1], "mean", np.mean(yr2009d790), "std", np.std(yr2009d790))
print(2010, labelsc[1], "mean", np.mean(yr2010d790), "std", np.std(yr2010d790))
print(2011, labelsc[1], "mean", np.mean(yr2011d790), "std", np.std(yr2011d790))
print(2012, labelsc[1], " mean ", np.mean(yr2012d790), " std ", np.std(yr2012d790), " num patients low ", yr2012drg790lo, " num patients high ", yr2012drg790hi)

print(2008, labelsc[2], "mean", np.mean(yr2008d791), "std", np.std(yr2008d791))
print(2009, labelsc[2], "mean", np.mean(yr2009d791), "std", np.std(yr2009d791))
print(2010, labelsc[2], "mean", np.mean(yr2010d791), "std", np.std(yr2010d791))
print(2011, labelsc[2], "mean", np.mean(yr2011d791), "std", np.std(yr2011d791))
print(2012, labelsc[2], " mean ", np.mean(yr2012d791), " std ", np.std(yr2012d791), " num patients low ", yr2012drg791lo, " num patients high ", yr2012drg791hi)

print(2008, labelsc[3], "mean", np.mean(yr2008d792), "std", np.std(yr2008d792))
print(2009, labelsc[3], "mean", np.mean(yr2009d792), "std", np.std(yr2009d792))
print(2010, labelsc[3], "mean", np.mean(yr2010d792), "std", np.std(yr2010d792))
print(2011, labelsc[3], "mean", np.mean(yr2011d792), "std", np.std(yr2011d792))
print(2012, labelsc[3], " mean ", np.mean(yr2012d792), " std ", np.std(yr2012d792), " num patients low ", yr2012drg792lo, " num patients high ", yr2012drg792hi)

print(2008, labelsc[4], "mean", np.mean(yr2008d793), "std", np.std(yr2008d793))
print(2009, labelsc[4], "mean", np.mean(yr2009d793), "std", np.std(yr2009d793))
print(2010, labelsc[4], "mean", np.mean(yr2010d793), "std", np.std(yr2010d793))
print(2011, labelsc[4], "mean", np.mean(yr2011d793), "std", np.std(yr2011d793))
print(2012, labelsc[4], " mean ", np.mean(yr2012d793), " std ", np.std(yr2012d793), " num patients low ", yr2012drg793lo, " num patients high ", yr2012drg793hi)

print(2008, labelsc[5], "mean", np.mean(yr2008d794), "std", np.std(yr2008d794))
print(2009, labelsc[5], "mean", np.mean(yr2009d794), "std", np.std(yr2009d794))
print(2010, labelsc[5], "mean", np.mean(yr2010d794), "std", np.std(yr2010d794))
print(2011, labelsc[5], "mean", np.mean(yr2011d794), "std", np.std(yr2011d794))
print(2012, labelsc[5], " mean ", np.mean(yr2012d794), " std ", np.std(yr2012d794), " num patients low ", yr2012drg794lo, " num patients high ", yr2012drg794hi)

print(2008, labelsc[6], "mean", np.mean(yr2008d795), "std", np.std(yr2008d795))
print(2009, labelsc[6], "mean", np.mean(yr2009d795), "std", np.std(yr2009d795))
print(2010, labelsc[6], "mean", np.mean(yr2010d795), "std", np.std(yr2010d795))
print(2011, labelsc[6], "mean", np.mean(yr2011d795), "std", np.std(yr2011d795))
print(2012, labelsc[6], " mean ", np.mean(yr2012d795), " std ", np.std(yr2012d795), " num patients low ", yr2012drg795lo, " num patients high ", yr2012drg795hi)


# Prints the results of the lower loop

labels = ["370/765: C Section w/ CC/MCC", "371/766: C Section w/o CC/MCC", "374/767: Vag Del w/ Ster.", "375/768: Vag Del w/ OR Proc", "372/774: Vag Del w/ Comp Diag", "373/775: Vag Del w/o Comp Diag"]


mavg2008 = [0,0,0,0,0,0]
mavg2009 = [0,0,0,0,0,0]
mavg2010 = [0,0,0,0,0,0]
mavg2011 = [0,0,0,0,0,0]
mavg2012 = [0,0,0,0,0,0]

for i in range(0,6) :
    mavg2008[i] = mchrg2008[i]/moms2008[i]
    mavg2009[i] = mchrg2009[i]/moms2009[i]
    mavg2010[i] = mchrg2010[i]/moms2010[i]
    mavg2011[i] = mchrg2011[i]/moms2011[i]
    mavg2012[i] = mchrg2012[i]/moms2012[i]

for i in range(0,6):
    print(2008, labels[i], mavg2008[i])
    print(2009, labels[i], mavg2009[i])
    print(2010, labels[i], mavg2010[i])
    print(2011, labels[i], mavg2011[i])
    print(2012, labels[i], mavg2012[i])


# Prints the variances:

print(2008, labels[0], "mean", np.mean(myr2008d765), "std", np.std(myr2008d765))
print(2009, labels[0], "mean", np.mean(myr2009d765), "std", np.std(myr2009d765))
print(2010, labels[0], "mean", np.mean(myr2010d765), "std", np.std(myr2010d765))
print(2011, labels[0], "mean", np.mean(myr2011d765), "std", np.std(myr2011d765))
print(2012, labels[0], "mean", np.mean(myr2012d765), "std", np.std(myr2012d765))

print(2008, labels[1], "mean", np.mean(myr2008d766), "std", np.std(myr2008d766))
print(2009, labels[1], "mean", np.mean(myr2009d766), "std", np.std(myr2009d766))
print(2010, labels[1], "mean", np.mean(myr2010d766), "std", np.std(myr2010d766))
print(2011, labels[1], "mean", np.mean(myr2011d766), "std", np.std(myr2011d766))
print(2012, labels[1], "mean", np.mean(myr2012d766), "std", np.std(myr2012d766))

print(2008, labels[2], "mean", np.mean(myr2008d767), "std", np.std(myr2008d767))
print(2009, labels[2], "mean", np.mean(myr2009d767), "std", np.std(myr2009d767))
print(2010, labels[2], "mean", np.mean(myr2010d767), "std", np.std(myr2010d767))
print(2011, labels[2], "mean", np.mean(myr2011d767), "std", np.std(myr2011d767))
print(2012, labels[2], "mean", np.mean(myr2012d767), "std", np.std(myr2012d767))

print(2008, labels[3], "mean", np.mean(myr2008d768), "std", np.std(myr2008d768))
print(2009, labels[3], "mean", np.mean(myr2009d768), "std", np.std(myr2009d768))
print(2010, labels[3], "mean", np.mean(myr2010d768), "std", np.std(myr2010d768))
print(2011, labels[3], "mean", np.mean(myr2011d768), "std", np.std(myr2011d768))
print(2012, labels[3], "mean", np.mean(myr2012d768), "std", np.std(myr2012d768))

print(2008, labels[4], "mean", np.mean(myr2008d774), "std", np.std(myr2008d774))
print(2009, labels[4], "mean", np.mean(myr2009d774), "std", np.std(myr2009d774))
print(2010, labels[4], "mean", np.mean(myr2010d774), "std", np.std(myr2010d774))
print(2011, labels[4], "mean", np.mean(myr2011d774), "std", np.std(myr2011d774))
print(2012, labels[4], "mean", np.mean(myr2012d774), "std", np.std(myr2012d774))

print(2008, labels[5], "mean", np.mean(myr2008d775), "std", np.std(myr2008d775))
print(2009, labels[5], "mean", np.mean(myr2009d775), "std", np.std(myr2009d775))
print(2010, labels[5], "mean", np.mean(myr2010d775), "std", np.std(myr2010d775))
print(2011, labels[5], "mean", np.mean(myr2011d775), "std", np.std(myr2011d775))
print(2012, labels[5], "mean", np.mean(myr2012d775), "std", np.std(myr2012d775))

# These should be combined somehow, but in what way?


#
# for row in drgdata:
#     if row[2][0:8] == "HOSPITAL":
#         print(row)
