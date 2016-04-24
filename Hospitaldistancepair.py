#Hospital - distance Pair


import numpy as np
import csv
import pickle

# This is the basic info about all the zip codes in the country, including population
# and latitude and longitude.
data = [line.strip().split() for line in open("Gaz_zcta_national.txt", 'rb')]

#found these indices by searching - but a few are missing.  Def 88063, 88220, 88430
ind1 = 25859
ind2 = 27793
TXzip = data[ind1:ind2+1]

# Indices to find missing zips
i=0
j=0
k=0


for n in range(len(data)):
	if data[n][0]=='88430':
		i = n
	elif data[n][0]=='88220':
		j = n
	elif data[n][0]=='88063':
		k = n



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
	hospdata[i][14] = eval(hospdata[i][14])
	hospdata[i][15] = eval(hospdata[i][15])


for i in range(len(TXzip)):
	for j in range(1,len(hospdata)):
		dist = dfunc(TXzip[i][7], TXzip[i][8], hospdata[j][14], hospdata[j][15])
		if dist < 25:
			TXzip[i].append(hospdata[i][0:16])

 #Now the format here is that each hospital record is followed by the distance to the centroid of that
 #zip code.

with open('zipcenter_distances_from_hosps.csv', 'wb') as fp:
	print 'saving'
	a = csv.writer(fp, delimiter=',')
	a.writerows(TXzip)

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
