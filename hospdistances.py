# Hospdistances.py

'''

Created: 04 04 16
Version: 04 07 16

This takes the set of all addresses for TX hospital-year observations, computes distances
and then tracks facilities nearby by year.

- hosplatlong.py should be run first (some 2016 or later version)
- Then some merging needs to be done in stata to produce TX Unique Lats Lons.csv
- this should be done with TX Hospital Sets.do
- Then run the current file
- It will produce, for every hospital:
    :: a set of only the other facilities within 50 miles
    :: a count of each of the Level 1, 2 or 3 facilities at 0-5, 5-15 and 15-25 miles
- The result is saved again in TX Unique Lats Lons.csv

'''


import csv
import pickle
import requests
import urllib
from lxml import etree
import numpy as np

# Keep track of the columns where certain pieces of data are recorded:
addr_add = 8
city_add = 4
intensive_add = 6
soloint_add = 13
lat_add = 14
lon_add = 15
fid_add = 0
year_add = 7

hospdata = []

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Unique Lats Lons.csv', 'r') as f:
    a = csv.reader(f, delimiter=',')
    for row in a:
        hospdata.append(row)

# Numbers imported as unicode strings - replace with eval:

for row in range(1,len(hospdata)):
    hospdata[row][intensive_add] = eval(hospdata[row][intensive_add])
    hospdata[row][soloint_add] = eval(hospdata[row][soloint_add])
    hospdata[row][fid_add] = eval(hospdata[row][fid_add])
    hospdata[row][year_add] = eval(hospdata[row][year_add])

# If necessary - to start over again using different distances, use these lines
distdata = []
for row in hospdata:
    distdata.append(row[0:16])

# If not, use the following:
# distdata = hospdata


def dfunc (w,x,y,z): #let these be (w,x) = (lat, lon) and (y,z) = (lat, lon)
    rad = 3961  #this is the Miles radius.  Kilometers =  6371
    conv = np.pi/180 #this should convert from decimal degrees to radians = pi/180
    w = w*conv
    x = x*conv
    y = y*conv
    z = z*conv
    d = 2*rad*np.arcsin( np.sqrt( np.square(np.sin( (w - y)/2)) + np.cos(w)*np.cos(y)*np.square(np.sin((x - z)/2 )) ))
    return d

bin1 = '0-5 Miles' # column 17
bin1_c1 = 0 # Level 1 - col 18
bin1_c2 = 0 # Level 2 - col 19
bin1_c3 = 0 # Level 3 - col 20
bin2 = '5-15 Miles' # column 21
bin2_c1 = 0 # Level 1 - col 22
bin2_c2 = 0 # Level 2 - col 23
bin2_c3 = 0 # Level 3 - col 24
bin3 = '15-25 Miles' # column 25
bin3_c1 = 0 # Level 1 - col 26
bin3_c2 = 0 # Level 2 - col 27
bin3_c3 = 0 # Level 3 - col 28


for row in range(1,len(distdata)):
    # Track level 1, 2, 3 at distance 0 - 5
    distdata[row].append(bin1)
    distdata[row].append(bin1_c1)
    distdata[row].append(bin1_c2)
    distdata[row].append(bin1_c3)
    # Track level 1, 2, 3 at distance 5 - 15
    distdata[row].append(bin2)
    distdata[row].append(bin2_c1)
    distdata[row].append(bin2_c2)
    distdata[row].append(bin2_c3)
    # Track level 1, 2, 3 at distance 15 - 25
    distdata[row].append(bin3)
    distdata[row].append(bin3_c1)
    distdata[row].append(bin3_c2)
    distdata[row].append(bin3_c3)
    # own latitude and longitude
    lat = eval(distdata[row][lat_add]) # 14
    lon = eval(distdata[row][lon_add]) # 15
    for other in range(1, len(distdata)):
        if distdata[row][year_add] == distdata[other][year_add]:
            oth_dist = dfunc(lat, lon, eval(distdata[other][lat_add]), eval(distdata[other][lon_add]) )
            if (oth_dist < 25) & (oth_dist > 0): # will append records of those facilities in less than 25 miles, but greater than 0 (i.e., not self)
                if (oth_dist <= 5) & (oth_dist > 0):
                    if not (distdata[other][0] in distdata[row]):
                        if (distdata[other][intensive_add] == 1) & (distdata[other][soloint_add] == 0):
                            distdata[row][19] += 1 #19
                        elif (distdata[other][intensive_add] == 0) & (distdata[other][soloint_add] == 1):
                            distdata[row][18] += 1 # 18
                        elif (distdata[other][intensive_add] == 0) & (distdata[other][soloint_add] == 0):
                            distdata[row][17] += 1 #17
                        else:
                            print("What the hell...")
                            print(distdata[other][0:28])
                elif (oth_dist > 5) & (oth_dist <= 15):
                    if not (distdata[other][0] in distdata[row]):
                        if (distdata[other][intensive_add] == 1) & (distdata[other][soloint_add] == 0):
                            distdata[row][23] += 1
                        elif (distdata[other][intensive_add] == 0) & (distdata[other][soloint_add] == 1):
                            distdata[row][22] += 1
                        elif (distdata[other][intensive_add] == 0) & (distdata[other][soloint_add] == 0):
                            distdata[row][21] += 1
                        else:
                            print("What the hell...")
                            print(distdata[other][0:28])
                elif (oth_dist > 15):
                    if not (hospdata[other][0] in hospdata[row]):
                        if (distdata[other][intensive_add] == 1) & (distdata[other][soloint_add] == 0):
                            distdata[row][27] += 1
                        elif (distdata[other][intensive_add] == 0) & (distdata[other][soloint_add] == 1):
                            distdata[row][26] += 1
                        elif (distdata[other][intensive_add] == 0) & (distdata[other][soloint_add] == 0):
                            distdata[row][25] += 1
                        else:
                            print("What the hell...")
                            print(distdata[other][0:28])
                distdata[row].append(distdata[other][0])
                distdata[row].append(oth_dist)




with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX Unique Lats Lons.csv', 'w') as f:
    print ('saving')
    a = csv.writer(f, delimiter=',')
    a.writerows(distdata)
