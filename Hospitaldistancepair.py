#Hospital - distance Pair

# Original Version: 2014
# Current version: 04 24 16


import numpy as np
import csv
import pickle

# This is the basic info about all the zip codes in the country, including population
# and latitude and longitude.
data = [line.strip().split() for line in open("/Users/austinbean/Desktop/dynhosp/Gaz_zcta_national.txt", 'rb')]

TXzip = []
#found these indices by searching - but a few are missing.  Def 88063, 88220, 88430
ind1 = 25859
ind2 = 27793
TXzip = data[ind1:ind2+1]

# Indices to find missing zips
i=0
j=0
k=0


for n in range(len(data)):
	if data[n][0]==b'88430':
		i = n
		print(i)
	elif data[n][0]==b'88220':
		j = n
		print(j)
	elif data[n][0]==b'88063':
		k = n
		print(k)



# This appends three missing zip codes starting with 8 to the data
TXzip.append(data[i])
TXzip.append(data[j])
TXzip.append(data[k])

for i in range(len(TXzip)):
	for k in range(len(TXzip[i])):
		TXzip[i][k] = eval(TXzip[i][k])

def dfunc (w,x,y,z): #let these be (w,x) = (lat, lon) and (y,z) = (lat, lon)
	rad = 3961  #this is the Miles radius.  Kilometers =  6371
	conv = np.pi/180 #this should convert from decimal degrees to radians = pi/180
	w = w*conv
	x = x*conv
	y = y*conv
	z = z*conv
	d = 2*rad*np.arcsin( np.sqrt( np.square(np.sin( (w - y)/2)) + np.cos(w)*np.cos(y)*np.square(np.sin((x - z)/2 )) ))
	# Gives answers in Miles
	# May want a correction for ill-conditioning of arcsin.
	# Should use something like arcsin(min(np.asarray([1]*len(y)),...)
	return d


hospdata = []


with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Unique Lats Lons.csv', 'r') as f:
    a = csv.reader(f, delimiter=',')
    for row in a:
        hospdata.append(row)

for i in range(1,len(hospdata)):
	hospdata[i][0] = eval(hospdata[i][0])
	hospdata[i][14] = eval(hospdata[i][14])
	hospdata[i][15] = eval(hospdata[i][15])

# Add field names:
# TXzip[0].append(hospdata[0][0:10] + hospdata[0][13:16])

indices = [ i for i in range(11)] + [i for i in range(13, 16)]
labels = [ hospdata[0][i] for i in indices]+['distance']
labels[12] = 'latitude'
labels[13] = 'longitude'

for i in range(len(TXzip)):
	temp = set() 		# Only want to append when it hasn't already been appended!
	TXzip[i].append(0)  # count the number of appended hospitals
	TXzip[i].append(0)  # note in the second loop below that this is a 50 miler.
	for j in range(1,len(hospdata)):
		dist = dfunc(TXzip[i][7], TXzip[i][8], hospdata[j][14], hospdata[j][15])
		if (dist < 25) and (hospdata[j][0] not in temp):
			TXzip[i][9] += 1
			for k in indices:
				TXzip[i].append(hospdata[j][k])
			TXzip[i].append(dist)	# don't forget the distance itself
			temp.add(hospdata[j][0])

# For those without any hospitals within 25 miles, do it again but at 35 miles:

for i in range(len(TXzip)):
	temp = set()
	if TXzip[i][9] <= 1:
		TXzip[i][10] = 1 # Set an indicator to record that a 50 mile radius was used.
		for j in range(1, len(hospdata)):
			dist = dfunc(TXzip[i][7], TXzip[i][8], hospdata[j][14], hospdata[j][15])
			if (dist < 50) and (hospdata[j][0] not in temp):
				TXzip[i][9] += 1
				for k in indices:
					TXzip[i].append(hospdata[j][k])
				TXzip[i].append(hospdata[j][k])
				temp.add(hospdata[j][0])


# This section will generate a list of variable labels for the STATA do file.
init = 0
for i in range(len(TXzip)):
	if TXzip[i][9] > init:
		init = TXzip[i][9] # computes the length of the longest zip record

col_labels = data[0] +['Hospital Count', '50 mile range'] + init*labels

for el in range(0,len(col_labels)):
	print(" rename v"+str(el+1) + " "+ str(col_labels[el])+str(((el-11) // 15)+1))


# This section will save the resulting file.

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Zip Distances from Hosps.csv', 'w') as fp:
	print('saving')
	a = csv.writer(fp, delimiter=',')
	a.writerows(TXzip)


'''
# Now let's create the list of (zip code, hospital THCIC_ID, hospital name, distance) quadruples (still called triples).
for j in range(len(TXzip)):
	for i in range(19, len(TXzip[1]),  10): # New Order below: zip number, ziplat, ziplong, pricepoint, hospname, hosplat, hosplong, distance
		trip = [TXzip[j][0], TXzip[j][7], TXzip[j][8], TXzip[j][i], TXzip[j][i+1], TXzip[j][i+7], TXzip[j][i+8], TXzip[j][i+9] ]
		triples.append(trip)


with open('ziptriples.csv', 'wb') as fp:
	print 'saving'
	a = csv.writer(fp, delimiter=',')
	a.writerows(triples)

# Next tasks - Two hospitals each with the name Memorial Hospital and Trinity Medical Center.  Find them and distinguish them.
# Indices for the two memorial hospitals are: 351, 352.  Just take the THCIC_ID from the file and add that to the name or something
# that will make it easy to find them.  That doesn't mean that the problem is done though...
'''
