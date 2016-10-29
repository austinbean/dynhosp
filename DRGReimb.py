

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
                yr2008d789.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[0] += eval(row[14].replace(',','').strip())
                chrg2009[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d789.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[0] += eval(row[14].replace(',','').strip())
                chrg2010[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d789.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[0] += eval(row[14].replace(',','').strip())
                chrg2011[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d789.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[0] += eval(row[14].replace(',','').strip())
                chrg2012[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d789.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
        elif eval(row[9]) == 790:
            if '2008' in row[8]:
                pats2008[1] += eval(row[14].replace(',','').strip())
                chrg2008[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d790.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[1] += eval(row[14].replace(',','').strip())
                chrg2009[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d790.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[1] += eval(row[14].replace(',','').strip())
                chrg2010[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d790.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[1] += eval(row[14].replace(',','').strip())
                chrg2011[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d790.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[1] += eval(row[14].replace(',','').strip())
                chrg2012[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d790.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
        elif eval(row[9]) == 791:
            if '2008' in row[8]:
                pats2008[2] += eval(row[14].replace(',','').strip())
                chrg2008[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d791.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[2] += eval(row[14].replace(',','').strip())
                chrg2009[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d791.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[2] += eval(row[14].replace(',','').strip())
                chrg2010[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d791.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[2] += eval(row[14].replace(',','').strip())
                chrg2011[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d791.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[2] += eval(row[14].replace(',','').strip())
                chrg2012[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d791.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
        elif eval(row[9]) == 792:
            if '2008' in row[8]:
                pats2008[3] += eval(row[14].replace(',','').strip())
                chrg2008[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d792.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[3] += eval(row[14].replace(',','').strip())
                chrg2009[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d792.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[3] += eval(row[14].replace(',','').strip())
                chrg2010[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d792.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[3] += eval(row[14].replace(',','').strip())
                chrg2011[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d792.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[3] += eval(row[14].replace(',','').strip())
                chrg2012[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d792.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
        elif eval(row[9]) == 793:
            if '2008' in row[8]:
                pats2008[4] += eval(row[14].replace(',','').strip())
                chrg2008[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d793.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[4] += eval(row[14].replace(',','').strip())
                chrg2009[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d793.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[4] += eval(row[14].replace(',','').strip())
                chrg2010[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d793.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[4] += eval(row[14].replace(',','').strip())
                chrg2011[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d793.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[4] += eval(row[14].replace(',','').strip())
                chrg2012[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d793.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
        elif eval(row[9]) == 794:
            if '2008' in row[8]:
                pats2008[5] += eval(row[14].replace(',','').strip())
                chrg2008[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d794.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[5] += eval(row[14].replace(',','').strip())
                chrg2009[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d794.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[5] += eval(row[14].replace(',','').strip())
                chrg2010[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d794.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[5] += eval(row[14].replace(',','').strip())
                chrg2011[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d794.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[5] += eval(row[14].replace(',','').strip())
                chrg2012[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d794.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
        elif eval(row[9]) == 795:
            if '2008' in row[8]:
                pats2008[6] += eval(row[14].replace(',','').strip())
                chrg2008[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2008d795.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2009' in row[8]:
                pats2009[6] += eval(row[14].replace(',','').strip())
                chrg2009[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2009d795.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2010' in row[8]:
                pats2010[6] += eval(row[14].replace(',','').strip())
                chrg2010[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2010d795.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2011' in row[8]:
                pats2011[6] += eval(row[14].replace(',','').strip())
                chrg2011[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2011d795.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))
            elif '2012' in row[8]:
                pats2012[6] += eval(row[14].replace(',','').strip())
                chrg2012[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
                yr2012d795.append(eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())/eval(row[14].replace(',','').strip()))



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

myr2008d789 = []
myr2009d789 = []
myr2010d789 = []
myr2011d789 = []
myr2012d789 = []

myr2008d790 = []
myr2009d790 = []
myr2010d790 = []
myr2011d790 = []
myr2012d790 = []

myr2008d791 = []
myr2009d791 = []
myr2010d791 = []
myr2011d791 = []
myr2012d791 = []

myr2008d792 = []
myr2009d792 = []
myr2010d792 = []
myr2011d792 = []
myr2012d792 = []

myr2008d793 = []
myr2009d793 = []
myr2010d793 = []
myr2011d793 = []
myr2012d793 = []

myr2008d794 = []
myr2009d794 = []
myr2010d794 = []
myr2011d794 = []
myr2012d794 = []

myr2008d795 = []
myr2009d795 = []
myr2010d795 = []
myr2011d795 = []
myr2012d795 = []


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
            elif '2009' in row[8]:
                moms2009[0] += eval(row[14].replace(',','').strip())
                mchrg2009[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2010' in row[8]:
                moms2010[0] += eval(row[14].replace(',','').strip())
                mchrg2010[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2011' in row[8]:
                moms2011[0] += eval(row[14].replace(',','').strip())
                mchrg2011[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2012' in row[8]:
                moms2012[0] += eval(row[14].replace(',','').strip())
                mchrg2012[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
        elif eval(row[9]) == 766:
            if '2008' in row[8]:
                moms2008[1] += eval(row[14].replace(',','').strip())
                mchrg2008[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2009' in row[8]:
                moms2009[1] += eval(row[14].replace(',','').strip())
                mchrg2009[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2010' in row[8]:
                moms2010[1] += eval(row[14].replace(',','').strip())
                mchrg2010[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2011' in row[8]:
                moms2011[1] += eval(row[14].replace(',','').strip())
                mchrg2011[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2012' in row[8]:
                moms2012[1] += eval(row[14].replace(',','').strip())
                mchrg2012[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
        elif eval(row[9]) == 767:
            if '2008' in row[8]:
                moms2008[2] += eval(row[14].replace(',','').strip())
                mchrg2008[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2009' in row[8]:
                moms2009[2] += eval(row[14].replace(',','').strip())
                mchrg2009[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2010' in row[8]:
                moms2010[2] += eval(row[14].replace(',','').strip())
                mchrg2010[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2011' in row[8]:
                moms2011[2] += eval(row[14].replace(',','').strip())
                mchrg2011[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2012' in row[8]:
                moms2012[2] += eval(row[14].replace(',','').strip())
                mchrg2012[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
        elif eval(row[9]) == 768:
            if '2008' in row[8]:
                moms2008[3] += eval(row[14].replace(',','').strip())
                mchrg2008[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2009' in row[8]:
                moms2009[3] += eval(row[14].replace(',','').strip())
                mchrg2009[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2010' in row[8]:
                moms2010[3] += eval(row[14].replace(',','').strip())
                mchrg2010[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2011' in row[8]:
                moms2011[3] += eval(row[14].replace(',','').strip())
                mchrg2011[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2012' in row[8]:
                moms2012[3] += eval(row[14].replace(',','').strip())
                mchrg2012[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
        elif eval(row[9]) == 774:
            if '2008' in row[8]:
                moms2008[4] += eval(row[14].replace(',','').strip())
                mchrg2008[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2009' in row[8]:
                moms2009[4] += eval(row[14].replace(',','').strip())
                mchrg2009[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2010' in row[8]:
                moms2010[4] += eval(row[14].replace(',','').strip())
                mchrg2010[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2011' in row[8]:
                moms2011[4] += eval(row[14].replace(',','').strip())
                mchrg2011[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2012' in row[8]:
                moms2012[4] += eval(row[14].replace(',','').strip())
                mchrg2012[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
        elif eval(row[9]) == 775:
            if '2008' in row[8]:
                moms2008[5] += eval(row[14].replace(',','').strip())
                mchrg2008[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2009' in row[8]:
                moms2009[5] += eval(row[14].replace(',','').strip())
                mchrg2009[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2010' in row[8]:
                moms2010[5] += eval(row[14].replace(',','').strip())
                mchrg2010[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2011' in row[8]:
                moms2011[5] += eval(row[14].replace(',','').strip())
                mchrg2011[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
            elif '2012' in row[8]:
                moms2012[5] += eval(row[14].replace(',','').strip())
                mchrg2012[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())


# Prints the results of the upper loop:

labelsc= ["385/789 Neonates died or transferred", "790 Extreme immaturity/respiratory distress", "791 Prematurity w/ Major Problems", "792 Prematurity w/out major problems", "793 Full term neonate w/ major problems", "794 Neonate w/ other significant problems", "795 Normal Newborn"]

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



# Prints the results of the lower loop

labels = ["765: C Section w/ CC/MCC", "766: C Section w/o CC/MCC", "767: Vag Del w/ Ster.", "768: Vag Del w/ OR Proc", "774: Vag Del w/ Comp Diag", "775: Vag Del w/o Comp Diag"]


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



# These should be combined somehow, but in what way?


#
# for row in drgdata:
#     if row[2][0:8] == "HOSPITAL":
#         print(row)
