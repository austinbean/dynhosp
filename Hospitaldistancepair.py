#Hospital - distance Pair

# Original Version: 2014
# Current version: 05 03 16


import numpy as np
import csv
import pickle

# This is a collection of zip codes from a commercial data source.

var_zipdata = []
with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/US Census Zip Code Tabulation Areas 2010/zip_code_database.csv', 'r') as f:
    a = csv.reader(f, delimiter=',')
    for row in a:
        if row[5] == 'TX':
            var_zipdata.append(eval(row[0]))

txset = set(var_zipdata)

# This is the basic info about all the zip codes in the country, including population
# and latitude and longitude.
data = [line.strip().split() for line in open("/Users/austinbean/Google Drive/Annual Surveys of Hospitals/US Census Zip Code Tabulation Areas 2010/2010 US Census Zip Code Tabulation Areas.txt", 'r')]

TXzip = []

for row in data[1:len(data)]:
    try:
        eval(row[0]) #zip codes which start with a zero give a SyntaxError
        if eval(row[0]) in txset:
            TXzip.append([eval(x) for x in row])
    except SyntaxError:
        print(row[0])
        # do nothing - don't care about zeros

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

print("Creating Hospital Data")

# TX Unique Lats Long.csv created by TX Hospital Sets.do
# And is produced by hospdistances.py - adds the addresses

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Unique Lats Lons.csv', 'r') as f:
    a = csv.reader(f, delimiter=',')
    for row in a:
        hospdata.append(row)

# Column labels from hospdistances.py:  Originally generated in TX Hospital Sets.do
fid_add = 0
facility_add = 1
county_add = 2
countynumber_add = 3
city_add = 4
mstat_add = 5
totalbeds_add = 6
intensive_add = 7
year_add = 8
deliveries_add = 9
nonicutransfersout_add = 10
nicutransfersin_add = 11
nicutransfersout_add = 12
nfp_add = 13
addr_add = 14
locozip_add = 15
ftephys_add = 16
fteresidents_add = 17
fteothertrainees_add = 18
fteregnurses_add = 19
ftelpn_add = 20
ftenap_add = 21
fteother_add = 22
ftetotal_add = 23
yearsbirths_add = 24
firstyear_add = 25
lastyear_add = 26
soloint_add = 27
lat_add = 28
lon_add = 29

for i in range(1,len(hospdata)):
    hospdata[i][fid_add] = eval(hospdata[i][fid_add])
    hospdata[i][lat_add] = eval(hospdata[i][lat_add])
    hospdata[i][lon_add] = eval(hospdata[i][lon_add])

# Add field names:
# TXzip[0].append(hospdata[0][0:10] + hospdata[0][13:16])

indices = [ i for i in range(14)] + [16,23] + [soloint_add, lat_add, lon_add]
labels = [ hospdata[0][i] for i in indices]+['distance']
labels[17] = 'latitude'
labels[18] = 'longitude'

print("Appending Hospitals Less than 25 miles distant")

for i in range(len(TXzip)):
    temp = set()         # Only want to append when it hasn't already been appended!
    TXzip[i].append(0)  # count the number of appended hospitals
    TXzip[i].append(0)  # note in the second loop below that this is a 50 miler.
    for j in range(1,len(hospdata)):
        dist = dfunc(TXzip[i][7], TXzip[i][8], hospdata[j][lat_add], hospdata[j][lon_add])
        if (dist < 25) and (dist not in temp): # and (hospdata[j][fid_add] not in temp):
            TXzip[i][9] += 1
            for k in indices:
                TXzip[i].append(hospdata[j][k])
            TXzip[i].append(dist)    # don't forget the distance itself
            temp.add(dist)  # temp.add(hospdata[j][fid_add])
'''
The above changed a little bit - there appear to be some hospitals which have
different names, different fids, but the same locations (to 16 decimal places
of measured distance from the zip centroid).  Perhaps the distance should be the
unique item.  This is definitely a problem.  The code to change the contents of
temp = set() here contains both variants, one commented.  Either do:
dist not in temp
OR
hospdata[j][fid_add] not in temp
to exclude by distance or by fid.
'''


# we need a list of pairs: hospitals and zips, with a distance from one to the other for EVERY possible zip/hospital combination
# Probably need to start over again:

print("Computing all Hospital Distances from Zip Centers")

alldistances = [ [i[0], i[7], i[8]] for i in TXzip]

for i in range(len(alldistances)):
    temp = set()
    for j in range(1, len(hospdata)):
        dist = dfunc(alldistances[i][1], alldistances[i][2], hospdata[j][lat_add], hospdata[j][lon_add])
        if hospdata[j][0] not in temp:  # this "temp" does not reflect the distance correction above - keep it this way, I think
            alldistances[i].append(hospdata[j][fid_add])
            alldistances[i].append(hospdata[j][facility_add]) # Note that the indices here do NOT follow those in "indices" the list, but those in hospdata[0]
            alldistances[i].append(hospdata[j][soloint_add]) # solo Intermediate
            alldistances[i].append(hospdata[j][intensive_add]) # neo intensive
            alldistances[i].append(dist)
            temp.add(hospdata[j][0])

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Zip All Hospital Distances.csv', 'w') as fp:
    print('saving')
    a = csv.writer(fp, delimiter=',')
    a.writerows(alldistances)


# Given the above, now we do:
# For those without any hospitals within 25 miles, do it again but at 50 miles.
# For those with too many hospitals, take the 10 closest.

print("Finding Extra Hospitals for zip codes with only 1")

for i in range(len(TXzip)):
    temp = set()
    if TXzip[i][9] <= 1:
        TXzip[i][10] = 1 # Set an indicator to record that a 50 mile radius was used.
        for j in range(1, len(hospdata)):
            dist = dfunc(TXzip[i][7], TXzip[i][8], hospdata[j][lat_add], hospdata[j][lon_add])
            if (dist < 50) and (dist not in temp): # (hospdata[j][fid_add] not in temp):
                TXzip[i][9] += 1
                for k in indices:
                    TXzip[i].append(hospdata[j][k])
                TXzip[i].append(dist)
                temp.add(dist) # temp.add(hospdata[j][fid_add])

max_hosp = 10 # If there are more than 10 hospitals, take the closest 10.

# Ways to list names, latitudes or longitudes.  Not that the order is not the same
# allnames = sorted([ TXzip[i][x] for x in range(12, len(TXzip[i]), 20) ]) # Note THESE WILL NOT SORT IN THE SAME ORDER
# alllats = sorted([ TXzip[i][x] for x in range(28, len(TXzip[i]), 20) ])
# alllons = sorted([ TXzip[i][x] for x in range(29, len(TXzip[i]), 20) ])

print("Keeping the closest 10 for those zips with too many")
reclength = 20 # this is the length of a hospital record as modified using the selections in "indices"

for i in range(len(TXzip)):
    if TXzip[i][9] > max_hosp:
#        TXzip[i][9] = max_hosp
        sorteddist = sorted([ TXzip[i][x] for x in range(30, len(TXzip[i]), reclength) ]) # all distances to that zip code, sorted ascending
        accept = sorteddist[0:max_hosp] # take the first "max_hosp" of them
        reject = sorteddist[max_hosp:len(sorteddist)] # get rid of these
        for el in reject: # This has to be reject - these are being REMOVED
            try:
                ind = TXzip[i].index(el)
                TXzip[i] = TXzip[i][0:ind-reclength+1]+ TXzip[i][ind+1:len(TXzip[i])]
            except ValueError:
                print("at row ", i)
                print("distance ", el)


# This section will generate a list of variable labels for the STATA do file.

print("STATA do file variable labels")

init = 10
for i in range(len(TXzip)):
    if TXzip[i][9] > init:
        init = TXzip[i][9] # computes the length of the longest zip record

col_labels = data[0] +['Hospital Count', '50 mile range'] + init*labels

for el in range(0,len(col_labels)):
    print(" rename v"+str(el+1) + " "+ str(col_labels[el])+str(((el-11) //reclength)+1))


# This section will save the resulting file.

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Zip Distances from Hosps.csv', 'w') as fp:
    print('saving')
    a = csv.writer(fp, delimiter=',')
    a.writerows(TXzip)


'''
# Now we also want to make sure each person's choice set includes the actual hospital chosen.

Generate one set of records for each hospital-zip code pair in the whole state.

Zip code, Hospital FID, NeoIntensive or Intermediate, year, location, distance, etc.

This requires adding the fid to the records in the Inpatient PUDF.
'''

# This next loop generates for each zip code the closest hospital to the zip.

print("Finding the closest hospital to each zip code")

for i in range(len(TXzip)):
    sorteddist = sorted([ TXzip[i][x] for x in range(30, len(TXzip[i]), reclength) ]) # all distances to that zip code, sorted ascending
    try:
        target = sorteddist[0]
        ind = TXzip[i].index(target)
        TXzip[i] = TXzip[i][0:11]+ TXzip[i][ind-reclength+1:ind+1]
    except IndexError:
        print("*******")
        print("No Hospitals within 50 miles")
        print(TXzip[i][0])
        print(sorteddist)

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Zip Closest Hosp.csv', 'w') as fp:
    print('saving')
    a = csv.writer(fp, delimiter=',')
    a.writerows(TXzip)
