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
sumnames = []

for resfile in os.listdir(resdir):
	if resfile.endswith(".pkl") and resfile.startswith("res_"):
		datasets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		filenames.append(resfile)
	elif resfile.endswith(".pkl") and resfile.startswith("sumres_"):
		sumsets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		sumnames.append(resfile)

# sort names alphabetically
ordereddata = collections.OrderedDict(sorted(datasets.items()))
orderedsum = collections.OrderedDict(sorted(sumsets.items()))

# Sanity check: check whether the number of tested instances instances differs
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

# Sanity check: check for fails, let the dev know
for res in sumnames:
	if not sumsets[res][Fail] == 0:	
		print "--------------------------------------------------------------------------------------------------"
		print "Warning: There were some failed runs in the tests. This might influence the significance of the"
		print "comparisons! Recommendation: Check for memlimits, aborts, fails etc. in the tested GCG versions."
		print "--------------------------------------------------------------------------------------------------"
		break

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
		ax.text(item.get_x()+item.get_width()/2., position, '%d' % int(height), ha='center', size='xx-small')
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
plt.title("Runtime per version")
plt.xlabel("GCG Version")
plt.ylabel("Runtime in seconds")
bars = plt.bar(range(len(runtime)), runtime.values(), align='center')
plt.xticks(range(len(runtime)), runtime.keys(), rotation=90)
setbarplotparams(highesttime)
plt.savefig(outdir + '/runtimes.pdf')				# name of image

# 3) Plot runtime comparison
# Calculate version-to-version speedup

items = list(runtime.items())
if len(items) < 2:
	print "Enter more than one GCG version to generate a runtime comparison plot."
else:
	highestdiff = 0
	runtimecomp = collections.OrderedDict()
	cumulative = collections.OrderedDict()
	for i in range(len(items)):
		if i > 0:
			# from the second item on calculate the version speed differences
			name = items[i-1][0] + '\n->\n' + items[i][0]
			diff = float(items[i-1][1]) + float(items[i][1])			
			runtimecomp[name] = diff
			if diff > highestdiff:
				highestdiff = diff
			# for the first one set initial cumulative value
			if i == 1:
				cumulative[name] = diff
			# for all following add the last value to current diff
			else:
				cumitems = list(cumulative.items())
				cumsum = cumitems[len(cumitems)-1][1]
				cumulative[name] = diff + cumsum

	# first plot version-to-version comparison bars
	fig, ax1 = plt.subplots()
	bar1 = ax1.bar(range(len(runtimecomp)), runtimecomp.values(), color='b')
	plt.xticks(range(len(runtimecomp)), runtimecomp.keys(), rotation=90)
	plt.tick_params(axis='x', which='major', labelsize=5)
	ax1.set_ylabel('Speedup in seconds', color='b')
	ax1.tick_params('y', colors='b')

	# then plot cumulative speedup
	ax2 = ax1.twinx()
	ax2.plot(range(len(runtimecomp)), cumulative.values(), 'r-')
	ax2.set_ylabel('Cumulative Speedup in seconds', color='r')
	ax2.tick_params('y', colors='r')
	
	fig.tight_layout()
	plt.title("Version-to-version runtime comparison")
	plt.savefig(outdir + '/runtimecomparison.pdf')			# name of image
