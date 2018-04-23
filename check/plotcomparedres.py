#!/usr/bin/env python
# Comparison plots for res files. 
# Takes res files in pickles format (aquired by parseres.py).
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import collections

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
	if not os.path.exists(outdir):
	    os.makedir(outdir)

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

# sort names alphabetically
ordereddata = collections.OrderedDict(sorted(datasets.items()))

# Check whether the number of tested instances instances differs (sanity check)
ninstances = -1
printwarning = False
for res in filenames:
	if ninstances == -1:
		ninstances = ordereddata[res].shape[0]	# count rows
	else:
		if ninstances != ordereddata[res].shape[0]:
			printwarning = True
if printwarning == True:
	print "--------------------------------------------------------------------------------------------------"
	print "Warning: Not all tests had the same number of instances."
	print "Did you enter more than one testset? Did all tested versions have access to all testset instances?"
	print "--------------------------------------------------------------------------------------------------"


# Get some statistics for each res file (first in temp dicts that will later be sorted)
maxstringlen = 12 # TODO make this number flexible
tempfails = {}
highestfails = 0
tempruntime = {}
highesttime = 0

for key in ordereddata.keys():
	# crop the filenames (keys in ordereddata) by removing res_ ... .pkl and add linebreak for very long keys
	croppedkey = key.split('/')[-1].replace('res_', '').replace('.pkl', '')
	if len(croppedkey) > maxstringlen:
		charlist = list(croppedkey)
		ninserts = len(charlist)/maxstringlen
		while ninserts > 0:
			charlist.insert(ninserts*maxstringlen, '\n')
			ninserts = ninserts - 1
		croppedkey = ''.join(charlist)
	# get amount of failed instances
	tempfails[croppedkey] = sumsets['sum' + key].loc['Fail']
	if tempfails[croppedkey] > highestfails:
		highestfails = tempfails[croppedkey]
	# get runtime
	tempruntime[croppedkey] = 0.0
	for time in ordereddata[key]['TotalTime']:
		tempruntime[croppedkey] = tempruntime[croppedkey] + float(time)
		if highesttime < tempruntime[croppedkey]:
			highesttime = tempruntime[croppedkey]

# order statistics by keys
fails = collections.OrderedDict(sorted(tempfails.items()))
runtime = collections.OrderedDict(sorted(tempruntime.items()))

# add a settings function
def setbarplotparams(highestbar):
	plt.ylim(ymin=0)
	if highestbar >= 10:
		valymax = highestbar+(highestbar/10)
	# guarantee that max y value is set to more than highest value
	elif highestbar == 0: 
		valymax = highestbar+2
	else:
		valymax = highestbar+1
	plt.ylim(ymax=valymax)	
	ax.grid(True,axis='y')
	plt.tight_layout()
	plt.tick_params(axis='x', which='major', labelsize=7)

	for item in bars:
		height = item.get_height()
		position = 1
		if highestbar > 0:
			position = height + (int(float(highestbar))/100)
		ax.text(item.get_x()+item.get_width()/2., position, '%ds' % int(height), ha='center', size='xx-small')
	return;

# 1) Plot how many instances were unsolved per version
fig = plt.figure()
ax = plt.axes()        
plt.title("Number of unsolved instances")
plt.xlabel("GCG Version")
bars = plt.bar(range(len(fails)), fails.values(), align='center')
plt.xticks(range(len(fails)), fails.keys(), rotation=90)
setbarplotparams(int(float(highestfails)))
plt.savefig(outdir + '/failcomparison.pdf')			# name of image

# 2) Plot runtime per version
fig = plt.figure()
ax = plt.axes()        
plt.title("Runtime comparison")
plt.xlabel("GCG Version")
plt.ylabel("Runtime in seconds")
bars = plt.bar(range(len(runtime)), runtime.values(), align='center')
plt.xticks(range(len(runtime)), runtime.keys(), rotation=90)
setbarplotparams(highesttime)
plt.savefig(outdir + '/runtimecomparison.pdf')			# name of image
