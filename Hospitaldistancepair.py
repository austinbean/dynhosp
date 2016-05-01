#Hospital - distance Pair

# Original Version: 2014
# Current version: 04 29 16


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

print("Appending Hospitals Less than 25 miles distant")

for i in range(len(TXzip)):
    temp = set()         # Only want to append when it hasn't already been appended!
    TXzip[i].append(0)  # count the number of appended hospitals
    TXzip[i].append(0)  # note in the second loop below that this is a 50 miler.
    for j in range(1,len(hospdata)):
        dist = dfunc(TXzip[i][7], TXzip[i][8], hospdata[j][14], hospdata[j][15])
        if (dist < 25) and (hospdata[j][0] not in temp):
            TXzip[i][9] += 1
            for k in indices:
                TXzip[i].append(hospdata[j][k])
            TXzip[i].append(dist)    # don't forget the distance itself
            temp.add(hospdata[j][0])

# we need a list of pairs: hospitals and zips, with a distance from one to the other for EVERY possible zip/hospital combination
# Probably need to start over again:

print("Computing all Hospital Distances from Zip Centers")

alldistances = [ [i[0], i[7], i[8]] for i in TXzip]

for i in range(len(alldistances)):
    temp = set()
    for j in range(1, len(hospdata)):
        dist = dfunc(alldistances[i][1], alldistances[i][2], hospdata[j][14], hospdata[j][15])
        if hospdata[j][0] not in temp:
            alldistances[i].append(hospdata[j][0])
            alldistances[i].append(hospdata[j][1]) # Note that the indices here do NOT follow those in "indices" the list, but those in hospdata[0]
            alldistances[i].append(hospdata[j][13]) # solo Intermediate
            alldistances[i].append(hospdata[j][6]) # neo intensive
            alldistances[i].append(dist)
            temp.add(hospdata[j][0])

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Zip All Hospital Distances.csv', 'w') as fp:
    print('saving')
    a = csv.writer(fp, delimiter=',')
    a.writerows(alldistances)



# For those without any hospitals within 25 miles, do it again but at 50 miles.
# For those with too many hospitals, take the 10 closest.

print("Finding Extra Hospitals for zip codes with only 1")

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
                TXzip[i].append(dist)
                temp.add(hospdata[j][0])

max_hosp = 10 # If there are more than 10 hospitals, take the closest 10.

print("Keeping the closest 10 for those zips with too many")

for i in range(len(TXzip)):
    if TXzip[i][9] > max_hosp:
#        TXzip[i][9] = max_hosp
        sorteddist = sorted([ TXzip[i][x] for x in range(25, len(TXzip[i]), 15) ]) # all distances to that zip code, sorted ascending
        accept = sorteddist[0:max_hosp] # take the first "max_hosp" of them
        reject = sorteddist[max_hosp:len(sorteddist)] # get rid of these
        for el in reject:
            try:
                ind = TXzip[i].index(el)
                TXzip[i] = TXzip[i][0:ind-14]+ TXzip[i][ind+1:len(TXzip[i])]
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
    print(" rename v"+str(el+1) + " "+ str(col_labels[el])+str(((el-11) // 15)+1))


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
    sorteddist = sorted([ TXzip[i][x] for x in range(25, len(TXzip[i]), 15) ]) # all distances to that zip code, sorted ascending
    try:
        target = sorteddist[0]
        ind = TXzip[i].index(target)
        TXzip[i] = TXzip[i][0:11]+ TXzip[i][ind-14:ind+1]
    except IndexError:
        print("*******")
        print("No Hospitals within 50 miles")
        print(TXzip[i][0])
        print(sorteddist)

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Zip Closest Hosp.csv', 'w') as fp:
    print('saving')
    a = csv.writer(fp, delimiter=',')
    a.writerows(TXzip)
