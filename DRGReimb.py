

import csv

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
				pats2008[0] += eval(row[14].replace(',','').strip())
				chrg2008[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2009' in row[8]:
				pats2009[0] += eval(row[14].replace(',','').strip())
				chrg2009[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2010' in row[8]:
				pats2010[0] += eval(row[14].replace(',','').strip())
				chrg2010[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2011' in row[8]:
				pats2011[0] += eval(row[14].replace(',','').strip())
				chrg2011[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2012' in row[8]:
				pats2012[0] += eval(row[14].replace(',','').strip())
				chrg2012[0] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
		elif eval(row[9]) == 790:
			if '2008' in row[8]:
				pats2008[1] += eval(row[14].replace(',','').strip())
				chrg2008[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2009' in row[8]:
				pats2009[1] += eval(row[14].replace(',','').strip())
				chrg2009[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2010' in row[8]:
				pats2010[1] += eval(row[14].replace(',','').strip())
				chrg2010[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2011' in row[8]:
				pats2011[1] += eval(row[14].replace(',','').strip())
				chrg2011[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2012' in row[8]:
				pats2012[1] += eval(row[14].replace(',','').strip())
				chrg2012[1] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
		elif eval(row[9]) == 791:
			if '2008' in row[8]:
				pats2008[2] += eval(row[14].replace(',','').strip())
				chrg2008[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2009' in row[8]:
				pats2009[2] += eval(row[14].replace(',','').strip())
				chrg2009[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2010' in row[8]:
				pats2010[2] += eval(row[14].replace(',','').strip())
				chrg2010[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2011' in row[8]:
				pats2011[2] += eval(row[14].replace(',','').strip())
				chrg2011[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2012' in row[8]:
				pats2012[2] += eval(row[14].replace(',','').strip())
				chrg2012[2] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
		elif eval(row[9]) == 792:
			if '2008' in row[8]:
				pats2008[3] += eval(row[14].replace(',','').strip())
				chrg2008[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2009' in row[8]:
				pats2009[3] += eval(row[14].replace(',','').strip())
				chrg2009[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2010' in row[8]:
				pats2010[3] += eval(row[14].replace(',','').strip())
				chrg2010[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2011' in row[8]:
				pats2011[3] += eval(row[14].replace(',','').strip())
				chrg2011[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2012' in row[8]:
				pats2012[3] += eval(row[14].replace(',','').strip())
				chrg2012[3] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
		elif eval(row[9]) == 793:
			if '2008' in row[8]:
				pats2008[4] += eval(row[14].replace(',','').strip())
				chrg2008[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2009' in row[8]:
				pats2009[4] += eval(row[14].replace(',','').strip())
				chrg2009[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2010' in row[8]:
				pats2010[4] += eval(row[14].replace(',','').strip())
				chrg2010[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2011' in row[8]:
				pats2011[4] += eval(row[14].replace(',','').strip())
				chrg2011[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2012' in row[8]:
				pats2012[4] += eval(row[14].replace(',','').strip())
				chrg2012[4] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
		elif eval(row[9]) == 794:
			if '2008' in row[8]:
				pats2008[5] += eval(row[14].replace(',','').strip())
				chrg2008[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2009' in row[8]:
				pats2009[5] += eval(row[14].replace(',','').strip())
				chrg2009[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2010' in row[8]:
				pats2010[5] += eval(row[14].replace(',','').strip())
				chrg2010[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2011' in row[8]:
				pats2011[5] += eval(row[14].replace(',','').strip())
				chrg2011[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2012' in row[8]:
				pats2012[5] += eval(row[14].replace(',','').strip())
				chrg2012[5] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
		elif eval(row[9]) == 795:
			if '2008' in row[8]:
				pats2008[6] += eval(row[14].replace(',','').strip())
				chrg2008[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2009' in row[8]:
				pats2009[6] += eval(row[14].replace(',','').strip())
				chrg2009[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2010' in row[8]:
				pats2010[6] += eval(row[14].replace(',','').strip())
				chrg2010[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2011' in row[8]:
				pats2011[6] += eval(row[14].replace(',','').strip())
				chrg2011[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())
			elif '2012' in row[8]:
				pats2012[6] += eval(row[14].replace(',','').strip())
				chrg2012[6] += eval(row[14].replace(',','').strip())*eval(row[15].replace('$','').replace(',','').strip())



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






for row in drgdata:
	if row[2][0:8] == "HOSPITAL":
		print(row)
