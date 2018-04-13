#!/usr/bin/env python
# This script reads *.res files created by GCG's make test and parses it into a pandas dataframe.
# Based on the parseout.py script
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# check command line arguments
if len(sys.argv) < 2:
	sys.exit("Usage: ./parseout.py RESFILE")

# array for all processed lines
linearray = [] 
columns = []

# checkout outfile
resfile = sys.argv[1]
fh = open(resfile, 'r')

# write solving information of testset into 
for line in fh:
	# get data column names (and add 'status' for the last column)
	if line.startswith("Name   "):
		line = " ".join(line.split())
		line = line.replace(" ", "")
		columns = line.split("|")
		columns.append("status")

	# get all data lines of first table
	elif '----------' not in line and line not in ['\n', '\r\n'] and not line.startswith("Name   ") and not line.startswith(" ") and not line.startswith("@"):
		line = " ".join(line.split())
		row = line.split(" ")
		# extra spaces should only be in the status column, join these
		if len(columns) < len(row):
			row[len(columns)-1:len(row)] = [''.join(row[len(columns)-1:len(row)])]
				
		linearray.append(row)

# there might be empty items in our coulumns list, remove these
index = 0
while index < len(columns):
	if columns[index] == "":
		del columns[index]
	else:
		index = index + 1

#store data into panda dataframe & save it as pickle
df = pd.DataFrame(columns=columns, data=linearray)
df.to_pickle('pickles/' + 'res_' + resfile.split('/')[-1].replace('.res', '.pkl'))
