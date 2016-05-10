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


for i in range(0, len(hospdata)):
	if hospdata[i][len(hospdata[i])-1] != 'FOUND':
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
			if (len(lat_set) >0) & (len(lon_set) > 0):
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
			time.sleep(0.55)
	elif hospdata[i][len(hospdata[i])-1] == 'FOUND':
		curr_add = hospdata[i][addr_add]

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
