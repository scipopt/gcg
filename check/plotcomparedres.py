#!/usr/bin/env python
# Comparison plots for res files. 
# Takes res files in pickles format (aquired by parseres.py).
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# check command line arguments
if len(sys.argv) < 2:
	sys.exit("Usage: ./plotcomparedres.py PKLDIR OUTPUTDIR (where OUTPUTDIR is optional)")

# get parameters
resdir = sys.argv[1]

outdirset = False
outdir = "pickles/"
if len(sys.argv) > 1:
	outdir = sys.argv[1]
	outdirset = True

# TODO put parsing here

# 1) Plot how many instances were unsolved per version

# Get premade res data
datasets = {}
sumsets = {}
filenames = []

for resfile in os.listdir(resdir):
	if resfile.endswith(".pkl") and resfile.startswith("res_"):
		datasets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		filenames.append(resfile)
	elif resfile.endswith(".pkl") and resfile.startswith("sumres_"):
		sumsets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))

# Check whether the number of tested instances instances differs (sanity check)
ninstances = -1
printwarning = False
for res in filenames:
	if ninstances == -1:
		ninstances = datasets[res].shape[0]	# count rows
	else:
		if ninstances != datasets[res].shape[0]:
			printwarning = True
if printwarning == True:
	print "--------------------------------------------------------------------------------------------------"
	print "Warning: Not all tests had the same number of instances."
	print "Did you enter more than one testset? Did all tested versions have access to all testset instances?"
	print "--------------------------------------------------------------------------------------------------"

# Count number of not "ok"/"solved"/"solved not verified" instances (failed/aborted/timeout/...) for each res file
fails = {}
highestfails = 0
maxstringlen = 12 # TODO make this number flexible

for key in datasets.keys():
	# crop the filenames (keys in datasets) by removing res_ ... .pkl and add linebreak for very long keys
	croppedkey = key.split('/')[-1].replace('res_', '').replace('.pkl', '')
	if len(croppedkey) > maxstringlen:
		charlist = list(croppedkey)
		ninserts = len(charlist)/maxstringlen
		while ninserts > 0:
			charlist.insert(ninserts*maxstringlen, '\n')
			ninserts = ninserts - 1
		croppedkey = ''.join(charlist)
	# get amount of failed instances
	fails[croppedkey] = sumsets['sum' + key].loc['Fail']
	if fails[croppedkey] > highestfails:
		highestfails = int(float(fails[croppedkey]))

# Plot results
fig = plt.figure()
ax = plt.axes()        
plt.title("Number of unsolved instances")
plt.xlabel("GCG Version")
plt.ylim(ymin=0)
if highestfails >= 10:
	valymax = highestfails+(highestfails/10)
# guarantee that max y value is set to more than highest value
elif highestfails == 0: 
	valymax = highestfails+2
else:
	valymax = highestfails+1
plt.ylim(ymax=valymax)	
ax.grid(True,axis='y')
bars = plt.bar(range(len(fails)), fails.values(), align='center')
plt.xticks(range(len(fails)), fails.keys(), rotation=90)
plt.tight_layout()
plt.tick_params(axis='x', which='major', labelsize=7)

for item in bars:
        height = item.get_height()
	if height == 0:
		ax.text(item.get_x()+item.get_width()/2., 1, '%d' % int(height), ha='center')
	else:
        	ax.text(item.get_x()+item.get_width()/2., 1.01*height, '%d' % int(height), ha='center')

plt.savefig(outdir + '/failcomparison.pdf')			# name of image
