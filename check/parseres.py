#!/usr/bin/env python3
# This script reads *.res files created by GCG's make test and parses it into a pandas dataframe.
# Based on the parseout.py script
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# check command line arguments
if len(sys.argv) < 2:
	sys.exit("Usage: ./parseout.py RESFILE OUTPUTDIR (where OUTPUTDIR is optional)")

# array for all processed lines of data table
linearray = [] 
columns = []

# array for (first) line of summary table
sumline = []
sumcolumns = []

# variable for timelimit
timelimit = 0

# get parameters
resfile = sys.argv[1]

outdirset = False
if len(sys.argv) > 2:
	outdir = sys.argv[2]
	outdirset = True
	if not os.path.exists(outdir):
	    os.makedirs(outdir)

# checkout outfile
fh = open(resfile, 'r')

# write solving information of testset, differentiate between big result table and small summary table
isdatatable = True
for line in fh:
	# get data column names (and add 'status' for the last column)
	if isdatatable and line.startswith("Name   "):
		line = " ".join(line.split())
		line = line.replace(" ", "")
		columns = line.split("|")
		columns.append("status")

	# get all data lines of first table
	elif isdatatable and '----------' not in line and line not in ['\n', '\r\n'] and not line.startswith("Name   ") and not line.startswith(" ") and not line.startswith("@"):
		line = " ".join(line.split())
		row = line.split(" ")
		linearray.append(row)

	# when summary table starts get column names and set isdatatable to False
	elif line.startswith("  Cnt "):
		isdatatable = False
		line = " ".join(line.split())
		sumcolumns = line.split(" ")
	
	# if isdatatable is False and the summary values line is reached get the data line
	elif not isdatatable and line.startswith("  "):
		line = " ".join(line.split())
		sumline = line.split(" ")
	# if the position is beyond all the other cases get the timelimit & finish reading
	elif not isdatatable and line.startswith("@02 timelimit: "):
		prefix, strtimelimit, timelimit = line.split(" ")
		timelimit = timelimit.replace("\n", "")
		break


# there might be empty items in our columns list, remove these
index = 0
while index < len(columns):
	if columns[index] == "":
		del columns[index]
	else:
		index = index + 1

# rename the last time column to TotalTime
for i, label in reversed(list(enumerate(columns))):
	if label == 'Time':
		columns[i] = 'TotalTime'
		break

# as the status column might contain spaces join all additional items
for row in linearray:
	if len(row) > len(columns):
		row[len(columns)-1:len(row)] = [''.join(row[len(columns)-1:len(row)])]

# store data into panda dataframe & save it as pickle
df = pd.DataFrame(columns=columns, data=linearray)
sumdf = pd.Series(index=sumcolumns, data=sumline)
timelimitdata = {'timelimit': [timelimit]}
timelimitdf = pd.DataFrame(data=timelimitdata)

if outdirset:
	df.to_pickle(outdir + '/' + 'res_' + resfile.split('/')[-1].replace('.res', '.pkl'))
	sumdf.to_pickle(outdir + '/' + 'sumres_' + resfile.split('/')[-1].replace('.res', '.pkl'))
	timelimitdf.to_pickle(outdir + '/' + 'timelimit_' + resfile.split('/')[-1].replace('.res', '.pkl'))
else:
	df.to_pickle('pickles/' + 'res_' + resfile.split('/')[-1].replace('.res', '.pkl'))
	sumdf.to_pickle('pickles/' + 'sumres_' + resfile.split('/')[-1].replace('.res', '.pkl'))
	timelimitdf.to_pickle('pickles/' + 'timelimit_' + resfile.split('/')[-1].replace('.res', '.pkl'))
