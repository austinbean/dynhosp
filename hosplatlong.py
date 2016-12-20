#hosplatlong.py

'''
Version 05 03 16 -
This now uses time.sleep() to avoid exceeding request rates.

- Takes TX All Addresses.csv as input.
- That file produced in TX Hospital Sets.do
- Outputs an updated version of the same, with addresses added.

'''

import csv
import pickle
import requests # This sometimes does not load.
import urllib
from lxml import etree
import numpy as np
import time

# Keep track of the columns where certain pieces of data are recorded:
addr_add = 8
city_add = 4
intensive_add = 6
soloint_add = 13
lat_add = 14
lon_add = 15
fid_add = 0

hospdata = []

with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX All Addresses.csv', 'r') as f:
	a = csv.reader(f, delimiter=',')
	for row in a:
		hospdata.append(row)

for row in range(1,len(hospdata)):
	hospdata[row][intensive_add] = eval(hospdata[row][intensive_add])
	hospdata[row][soloint_add] = eval(hospdata[row][soloint_add])
	hospdata[row][fid_add] = eval(hospdata[row][fid_add])

url_start = 'https://maps.googleapis.com/maps/api/geocode/xml?'

urls = []
curr_add = " "

# To reset and start again: 
for i in range(0, len(hospdata)):
	if hospdata[i][len(hospdata[i])-1] == 'FOUND':
		hospdata[i][len(hospdata[i])-1] = ''

# find this fucking Cass County hospital.
for i in range(0,len(hospdata)):
	if hospdata[i][2] == "CASS":
		print(hospdata[i])



for i in range(0, len(hospdata)):
	if hospdata[i][9] == '': # missing zip code
		print("***************")
		print("FOUND MISSING", i, hospdata[i][0], hospdata[i][9])
		if (hospdata[i-1][0] == hospdata[i][0])&(hospdata[i-1][9]!=''):
			print("Reassign to previous")
			print(hospdata[i][0], hospdata[i-1][0], hospdata[i][9], hospdata[i-1][9])
			hospdata[i][9] = hospdata[i-1][9]
		elif (hospdata[i+1][0] == hospdata[i][0])&(hospdata[i+1][9]!=''):
			print("Reassign to following:")
			print(hospdata[i][0], hospdata[i+1][0], hospdata[i][9], hospdata[i+1][9])
			hospdata[i][9] = hospdata[i+1][9]
		else:
			print("Double Missing?  At row ", i)


misscount = 0
for i in range(0, len(hospdata)):
	if (hospdata[i][9] == '')|(hospdata[i][addr_add] == ''):
		misscount += 1
		print("****************")
		print("PRIOR")
		print(hospdata[i-1])
		print("ACTUAL")
		print(hospdata[i])
		print("SUBSEQUENT")
		print(hospdata[i+1])
		print("****************")
print("TOTAL MISSING", misscount)
			

for i in range(0, len(hospdata)):
	street = '' #reset these strings to be empty to avoid looking up in the wrong city.
	town = ''
	lat = ''
	lon = ''
	url_req = ''
	if hospdata[i][len(hospdata[i])-1] != 'FOUND':
		if not ((hospdata[i][addr_add] == '')): # search if the address is not missing?
			if curr_add == hospdata[i][addr_add]: 
				print(hospdata[i][addr_add])
				hospdata[i].append(lat)
				hospdata[i].append(lon)
				hospdata[i].append('FOUND')
				print( "doing nothing") 
			elif curr_add != hospdata[i][addr_add]:
				street = 'address='+urllib.parse.quote_plus(hospdata[i][addr_add])+','
				town = '+'+urllib.parse.quote_plus(hospdata[i][city_add])+','
				state = '+TX'
				url_req = url_start+street+town+state+'&sensor=false'
				urls.append(url_req)
				a = requests.get(url_req)
				page_xml = etree.XML(a.content)
				lat_set = page_xml.xpath('//location/lat')
				lon_set = page_xml.xpath('//location/lng')
				if (len(lat_set) >0) and (len(lon_set) > 0):
					lat = lat_set[0].text
					lon = lon_set[0].text
					print( hospdata[i][0], hospdata[i][1], lat, lon)
					# Saves lat and long to data file
					hospdata[i].append(lat)
					hospdata[i].append(lon)
					hospdata[i].append('FOUND')
					curr_add = hospdata[i][addr_add]
				else:
					curr_add = hospdata[i][addr_add]
				time.sleep(0.8)
		else:
			print("*******")
			print("Stuff is missing!")
			print(i)
	elif hospdata[i][len(hospdata[i])-1] == 'FOUND':
		curr_add = hospdata[i][addr_add]




# Ok - some of the URLS give completely different answers.  
# To check them, compute distance among constant fids I think.
for tx in urls:
	if "Atlanta" in tx:
		print(tx)
	elif "ATLANTA" in tx:
		print(tx)


def dfunc (w,x,y,z): #let these be (w,x) = (lat, lon) and (y,z) = (lat, lon)
    rad = 3961  #this is the Miles radius.  Kilometers =  6371
    conv = np.pi/180 #this should convert from decimal degrees to radians = pi/180
    w = w*conv
    x = x*conv
    y = y*conv
    z = z*conv
    d = 2*rad*np.arcsin( np.sqrt( np.square(np.sin( (w - y)/2)) + np.cos(w)*np.cos(y)*np.square(np.sin((x - z)/2 )) ))
    return d


dcnt = 0
for row in range(1,len(hospdata)):
	for row2 in range(1,len(hospdata)):
		if not(hospdata[row][8]==''):
			if not( hospdata[row2][8]==''):
				if hospdata[row][0] == hospdata[row2][0]:
					if dfunc(eval(hospdata[row][14]), eval(hospdata[row][15]), eval(hospdata[row2][14]), eval(hospdata[row2][15])) > 15:
						print("*******************")
						print("DISTANCE IS: ", dfunc(eval(hospdata[row][14]), eval(hospdata[row][15]), eval(hospdata[row2][14]), eval(hospdata[row2][15])))
						print(hospdata[row])
						print(hospdata[row2])
						dcnt += 1
print(dcnt)


# OK - there are several to fix here.  A couple where the distance is 10 and then 2 where it is more than 100 miles.  Fix those by hand.  

for row in range(1,len(hospdata)):
	if hospdata[row][0] == 670132:
		hospdata[row][14] = '33.1050940'
		hospdata[row][15] = '-94.1654897'
	elif hospdata[row][0] == 4813735:
		hospdata[row][14] = '29.329170'
		hospdata[row][15] = '-96.119180'


# Check who was not found:

for i in range(0, len(hospdata)):
	if hospdata[i][len(hospdata[i])-1] != 'FOUND':
		print("********************")
		print(i)
		print(hospdata[i])



with open('/Users/austinbean/Google Drive/Annual Surveys of Hospitals/TX All Addresses.csv', 'w') as f:
	print ('saving')
	a = csv.writer(f, delimiter=',')
	a.writerows(hospdata)




# just_locs = []
# for i in range(len(hospdata)):
# 	vals = [ hospdata[i][0], hospdata[i][1], hospdata[i][15], hospdata[i][19], hospdata[i][20]]
# 	just_locs.append(vals)
'''
output_file_2 = open('hosp_data_w_latlong.pkl', 'wb')

pickle.dump(hospdata, output_file_2)

with open('hosp_name_and_loc.csv', 'w') as fp:
	print 'saving'
	a = csv.writer(fp, delimiter=',')
	a.writerows(just_locs)
'''


# KEEP THIS VALUE - this is the google maps API key.  Might be useful
# &key={AIzaSyDLybqdd-5MrOrDJGN8D76paBJC7UJS9KY}
