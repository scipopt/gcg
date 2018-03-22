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

#  variables to perform line split with
columns = ['name','type','origconss','origvars','preconss','prevars','detector','decblocks','decrel','decmconss','decmvars','dualbound','primbound','gap','prcalls','prvars','prtime','mlptime','miters','node','time','status']

# array for all processed lines
linearray = [] 

# checkout outfile
resfile = sys.argv[1]
fh = open(resfile, 'r')

# write solving information of testset into 
for line in fh:
	# get all data lines of first table
	if '----------' not in line and line not in ['\n', '\r\n'] and not line.startswith("Name   ") and not line.startswith(" ") and not line.startswith("@"):
		line = " ".join(line.split())
		linearray.append(line.split(" "))

#store data into panda dataframe & save it as pickle
df = pd.DataFrame(columns=columns, data=linearray)
df.to_pickle('pickles/' + 'res_' + resfile.split('/')[-1].replace('.res', '.pkl'))
