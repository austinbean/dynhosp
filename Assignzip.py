#Assignzip.py

# Original 2014
# Version 04 23 16

import numpy as np
import csv
import pickle

# This is the basic info about all the zip codes in the country, including population
# and latitude and longitude.
data = [line.strip().split() for line in open("/Users/austinbean/Desktop/dynhosp/Gaz_zcta_national.txt", 'rb')]

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
	elif data[n][0]==b'88220':
		j = n
	elif data[n][0]==b'88063':
		k = n

# This appends three missing zip codes starting with 8 to the data
TXzip.append(data[i])
TXzip.append(data[j])
TXzip.append(data[k])

# These are just Texas zip codes identified as a continuous range from the file
# Gaz_zcta_national.  This is simply to test whether we have all of them using the next
# set of information from the DEC_10 file
set1 = set()
for i in TXzip:
    set1.add(i[0])


# This is a list from the census bureau of zip codes in Texas
f = open( '/Users/austinbean/Google Drive/Current Projects/Second year Paper/TXHospitalData/DEC_10_SF1_P1_with_ann.csv', 'rU')
data2 = list(tuple(rec) for rec in csv.reader(f, delimiter=','))
set2 = set()

# This loads those zips into a set
for j in data2:
	set2.add(j[1])

# Elements not in one or the other -
b = set2.difference(set1)
c = set1.difference(set2)
# Given the correction above the first three of these will be in TXzip where they belong
# In b:
# 88430 - NM zip code with a little strip in TX
# 88220 - NM zip code with a piece in TX
# 88063 - NM zip code with a small piece of El Paso
# 73949 - OK zip without any part of TX.  DON'T worry about this one.  Not a TX zip.


# This saves a pkl of all the TX zip codes, with
f = open('txzip.pkl', 'wb')
pickle.dump(TXzip, f)
f.close()


# Now with all of the TX zip codes in TXzip, just need to compute the distances from hospitals to those zips.
# This is nearly done w/ what has already been written in distances.py.
# But what do I want, in the end?
# (i.) For each hospital, a vector of distances to every zip code.
# (ii.) For each (zip code, hospital) pair, a distance
# (iii.) For all hospitals, a set of distances of one from the other.

# All TXzip entries are in unicode.  Evaluate as numbers:

for i in range(len(TXzip)):
	for k in range(len(TXzip[i])):
		TXzip[i][k] = eval(TXzip[i][k])

# Only TX zip list
justzips = []
for el in TXzip:
	justzips.append(el[0])

# Open hospital data array:

f_open = open('hosp_data_w_latlong.pkl', 'rb')
dist = pickle.load(f_open)

#NB: 655 facilities each with 20 entries in a row.


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


# For each hospital record, append distances to every zip plus lists of zips
for i in range(len(dist)):
	hosp_lat = eval(dist[i][7]) #lat and long are stored as unicodes - evaluate as numbers
	hosp_lon = eval(dist[i][8])
	for j in range(len(TXzip)):
		num = dfunc(hosp_lat, hosp_lon, TXzip[j][7], TXzip[j][8])
		dist[i].append(TXzip[j]) # Appends result to hospital record
		dist[i].append(num)
# Now the format is that there is a hospital record, after the record comes a zip code entry,
# then the distance from the pop-weighted centroid of the zip to the hospital.

with open('hosp_distances_from_zipcenters.csv', 'wb') as fp:
	print 'saving'
	a = csv.writer(fp, delimiter=',')
	a.writerows(dist)

with open('/Users/austinbean/Desktop/dynhosp/TXzipsonly.csv', 'w', newline='') as fw:
	aout = csv.writer(fw, delimiter=',')
	aout.writerow(justzips)












