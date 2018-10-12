#!/usr/bin/env python
# Comparison plots for res files. 
# Takes res files in pickles format (aquired by parseres.py).
# Python 2.7

import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import collections
import math

# check command line arguments
if len(sys.argv) < 2:
	sys.exit('Usage: ./plotcomparedres.py PKLDIR OUTPUTDIR (where OUTPUTDIR is optional)')

# -------------------------------------------------------------------------------------------------------------------------
# Get all necessary parameters and statistics
# -------------------------------------------------------------------------------------------------------------------------

# get parameters
resdir = sys.argv[1]

outdirset = False
outdir = 'pickles/'
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
	if resfile.endswith('.pkl') and resfile.startswith('res_'):
		datasets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		filenames.append(resfile)
	elif resfile.endswith('.pkl') and resfile.startswith('sumres_'):
		sumsets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		sumnames.append(resfile)

nversions = len(datasets)

# sort names alphabetically
ordereddata = collections.OrderedDict(sorted(datasets.items()))
orderedsum = collections.OrderedDict(sorted(sumsets.items()))

# Sanity check: check whether the number of tested instances differs
ninstances = -1
printwarning = False
for res in filenames:
	if ninstances == -1:
		ninstances = ordereddata[res].shape[0]	# count rows
	else:
		if ninstances != ordereddata[res].shape[0]:
			printwarning = True
if printwarning == True:
	print '--------------------------------------------------------------------------------------------------'
	print 'Warning: Not all tests had the same number of instances.'
	print 'Did you enter more than one testset? Did all tested versions have access to all testset instances?'
	print '--------------------------------------------------------------------------------------------------'
	ninstances = -1

# Sanity check: check for fails, let the dev know
for res in sumnames:
	if not sumsets[res]['Fail'] == '0':	
		print '--------------------------------------------------------------------------------------------------'
		print 'Warning: There were some failed runs in the tests. This might influence the significance of the'
		print 'comparisons! Recommendation: Check for memlimits, aborts, fails etc. in the tested GCG versions.'
		print '--------------------------------------------------------------------------------------------------'
		break

# -------------------------------------------------------------------------------------------------------------------------
# Get some statistics for each res file to be used in the plots
# -------------------------------------------------------------------------------------------------------------------------

maxstringlen = 12 # TODO make this number flexible

versions = []

fails = {}
aborts = {}
memlimits = {}
timeouts = {}
readerrors = {}

timefails = {}
timeaborts = {}
timememlimits = {}
timetimeouts = {}
timereaderrors = {}
timesolved = {}

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
	versions.append(croppedkey)
	# get fail types and their amounts
	fails[croppedkey] = 0
	aborts[croppedkey] = 0
	memlimits[croppedkey] = 0
	timeouts[croppedkey] = 0
	readerrors[croppedkey] = 0
	for status in ordereddata[key]['status']:
		if status == 'fail':
			fails[croppedkey] = fails[croppedkey] + 1
		elif status == 'abort':
			aborts[croppedkey] = aborts[croppedkey] + 1
		elif status == 'memlimit':
			memlimits[croppedkey] = memlimits[croppedkey] + 1
		elif status == 'timeout':
			timeouts[croppedkey] = timeouts[croppedkey] + 1 
		elif status == 'readerror':
			readerrors[croppedkey] = readerrors[croppedkey] + 1
	# get amount of failed instances (including limits)
	failamount = sumsets['sum' + key].loc['Fail']
	if int(failamount) + timeouts[croppedkey] + memlimits[croppedkey] > highestfails:
		highestfails = int(failamount) + timeouts[croppedkey] + memlimits[croppedkey]
	# get runtime
	tempruntime[croppedkey] = 0.0
	for time in ordereddata[key]['TotalTime']:
		tempruntime[croppedkey] = tempruntime[croppedkey] + float(time)
		if highesttime < tempruntime[croppedkey]:
			highesttime = tempruntime[croppedkey]
	# get runtime per status
	timefails[croppedkey] = 0
	timeaborts[croppedkey] = 0
	timememlimits[croppedkey] = 0
	timetimeouts[croppedkey] = 0
	timereaderrors[croppedkey] = 0
	timesolved[croppedkey] = 0
	tablelength = len(ordereddata[key].index)
	for i in range(0, tablelength-1):
		if ordereddata[key]['status'][i] == 'fail':
			timefails[croppedkey] = timefails[croppedkey] + float(ordereddata[key]['TotalTime'][i])
		elif ordereddata[key]['status'][i] == 'abort':
			timeaborts[croppedkey] = timeaborts[croppedkey] + float(ordereddata[key]['TotalTime'][i])
		elif ordereddata[key]['status'][i] == 'memlimit':
			timememlimits[croppedkey] = timememlimits[croppedkey] + float(ordereddata[key]['TotalTime'][i])
		elif ordereddata[key]['status'][i] == 'timeout':
			timetimeouts[croppedkey] = timetimeouts[croppedkey] + float(ordereddata[key]['TotalTime'][i])
		elif ordereddata[key]['status'][i] == 'readerror':
			timereaderrors[croppedkey] = timereaderrors[croppedkey] + float(ordereddata[key]['TotalTime'][i])
		else:
			timesolved[croppedkey] = timesolved[croppedkey] + float(ordereddata[key]['TotalTime'][i])
	# round up runtime per status
	timefails[croppedkey] = math.ceil(timefails[croppedkey])
	timeaborts[croppedkey] = math.ceil(timeaborts[croppedkey])
	timememlimits[croppedkey] = math.ceil(timememlimits[croppedkey])
	timetimeouts[croppedkey] = math.ceil(timetimeouts[croppedkey])
	timereaderrors[croppedkey] = math.ceil(timereaderrors[croppedkey])
	timesolved[croppedkey] = math.ceil(timesolved[croppedkey])

# order statistics by keys
versions = sorted(versions)

fails = collections.OrderedDict(sorted(fails.items()))
aborts = collections.OrderedDict(sorted(aborts.items()))
memlimits = collections.OrderedDict(sorted(memlimits.items()))
timeouts = collections.OrderedDict(sorted(timeouts.items()))
readerrors = collections.OrderedDict(sorted(readerrors.items()))
runtime = collections.OrderedDict(sorted(tempruntime.items()))

timefails = collections.OrderedDict(sorted(timefails.items()))
timeaborts = collections.OrderedDict(sorted(timeaborts.items()))
timememlimits = collections.OrderedDict(sorted(timememlimits.items()))
timetimeouts = collections.OrderedDict(sorted(timetimeouts.items()))
timereaderrors = collections.OrderedDict(sorted(timereaderrors.items()))
timesolved = collections.OrderedDict(sorted(timesolved.items()))

# -------------------------------------------------------------------------------------------------------------------------
# Add functions for often used parts of the plots
# -------------------------------------------------------------------------------------------------------------------------

# function to label bars
def labelbars(bars, highestbar):
	for item in bars:
		height = item.get_height()
		position = 1
		if highestbar > 0:
			position = height + (int(float(highestbar))/100)
		ax.text(item.get_x()+item.get_width()/2., position, '%d' % int(height), ha='center', size='xx-small')
	return;

# settings function for axis limits, layout and bar tests
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
	return;

# -------------------------------------------------------------------------------------------------------------------------
# 1) Plot how many instances were unsolved per version
# -------------------------------------------------------------------------------------------------------------------------

fig = plt.figure()
ax = plt.axes()        
plt.title('Number of unsolved instances')
plt.xlabel('GCG Version')

faildata = {'fails': fails.values(), 'aborts': aborts.values(), 'memlimits': memlimits.values(), 
	'timeouts': timeouts.values(), 'readerrors': readerrors.values()}
failbars = pd.DataFrame(data=faildata)
failbars.plot(kind='bar', stacked=True)

# label the stacked bars
for (ind,row) in failbars.iterrows():
	cumval = 0
	for column in failbars:
		val = row.loc[column]
		if not val == 0:
			cumval = cumval + val
			plt.annotate( val, xy = (ind, cumval - .5), horizontalalignment='center', verticalalignment='top',
				fontsize=8 )

plt.xticks(range(len(fails)), fails.keys(), rotation=90)
setbarplotparams(int(float(highestfails)))

ax1 = plt.subplot(111)
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5, fancybox=False, prop={'size': 'small'}, framealpha=1.0)

# if the number of instances differs
stringninstances = 'unknown or differed'
if ninstances >= 0:
	stringninstances = str(ninstances)

plt.figtext(.01,.01,'The total number of instances in the test (per version) was ' + stringninstances + '.', 
	size='x-small')

plt.savefig(outdir + '/failcomparison.pdf')			# name of image

# -------------------------------------------------------------------------------------------------------------------------
# 2) Plot runtime per version
# -------------------------------------------------------------------------------------------------------------------------

fig = plt.figure()
ax = plt.axes()        
plt.title('Runtime per version')
plt.xlabel('GCG Version')
plt.ylabel('Runtime in seconds')

bars = plt.bar(range(len(runtime)), runtime.values(), align='center')
plt.xticks(range(len(runtime)), runtime.keys(), rotation=90)
setbarplotparams(highesttime)
labelbars(bars, highesttime)

plt.savefig(outdir + '/runtimes.pdf')				# name of image

# -------------------------------------------------------------------------------------------------------------------------
# 3) Plot runtime comparison
# -------------------------------------------------------------------------------------------------------------------------

# Calculate version-to-version speedup

items = list(runtime.items())
if len(items) < 2:
	print 'Enter more than one GCG version to generate a runtime comparison plot.'
else:	
	highestdiff = 0
	lowestdiff = 0
	highestcum = 0
	lowestcum = 0
	runtimecomp = collections.OrderedDict()
	cumulative = collections.OrderedDict()
	for i in range(len(items)):
		diff = 0
		if i > 0:
			# from the second item on calculate the version speed differences
			name = items[i-1][0] + '\n->\n' + items[i][0]
			diff = int(round( float(items[i-1][1]) - float(items[i][1]) ))
			maxtime = max(items[i-1][1], items[i][1])
			diff = float(diff) / maxtime
			runtimecomp[name] = diff
			if diff > highestdiff:
				highestdiff = diff
			if diff < lowestdiff:
				lowestdiff = diff
			# for the first one set initial cumulative value
			if i == 1:
				cumulative[name] = diff
				highestcum = diff
				lowestcum = diff
			# for all following add the last value to current diff
			else:
				cumitems = list(cumulative.items())
				cumsum = cumitems[len(cumitems)-1][1]
				cumulative[name] = diff + cumsum
				if cumulative[name] > highestcum:
					highestcum = cumulative[name]
				elif cumulative[name] < lowestcum:
					lowestcum = cumulative[name]

	#determine axis min/max (leave space for bar labels
	axmin = lowestdiff+0.1*lowestdiff
	if lowestcum < axmin:
		axmin = lowestcum
	if axmin > 0:
		axmin = 0

	axmax = highestdiff+0.1*highestdiff
	if highestcum > axmax:
		axmax = highestcum
	if axmax < 0:
		axmax = 0

	# first plot version-to-version comparison bars
	fig, ax1 = plt.subplots()
	bar1 = ax1.bar(range(len(runtimecomp)), runtimecomp.values(), color='b')
	plt.xticks(range(len(runtimecomp)), runtimecomp.keys(), rotation=90)
	plt.tick_params(axis='x', which='major', labelsize=5)
	ax1.set_ylabel('Speedup factor', color='b')
	ax1.tick_params('y', colors='b')

	# make space far bar labels
	ax1.set_ylim(ymin=axmin, ymax=axmax)
	longestbar = highestdiff	
	if abs(lowestdiff) > longestbar:
		longestbar = abs(lowestdiff)

	labelbars(bar1, longestbar)

	# plot cumulative speedup if there is more than one bar
	if len(items) > 2:
		ax2 = ax1.twinx()
		ax2.plot(range(len(runtimecomp)), cumulative.values(), 'r-')
		ax2.set_ylabel('Cumulative speedup factor', color='r')
		ax2.tick_params('y', colors='r')
		ax2.axhline(y=0, color='xkcd:orange')
		ax2.set_ylim(ymin=axmin, ymax=axmax)
	
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	plt.title('Version-to-version runtime comparison')

	plt.savefig(outdir + '/runtimecomparison.pdf')			# name of image

# -------------------------------------------------------------------------------------------------------------------------
# 4) Plot relative time per status category (fail categories and solved)
# -------------------------------------------------------------------------------------------------------------------------

fig = plt.figure()
ax = plt.axes()        
plt.title('Runtime per solving status')
plt.xlabel('GCG Version')
plt.ylabel('Runtime in seconds')

# get times per status relative to runtime
reltimefails = collections.OrderedDict()
reltimeaborts = collections.OrderedDict()
reltimememlimits = collections.OrderedDict()
reltimetimeouts = collections.OrderedDict()
reltimereaderrors = collections.OrderedDict()
reltimesolved = collections.OrderedDict()

for vers in versions:
	reltime = float(timefails[vers]) / runtime[vers]
	reltimefails.update({vers : reltime})

	reltime = float(timeaborts[vers]) / runtime[vers]
	reltimeaborts.update({vers : reltime})

	reltime = float(timememlimits[vers]) / runtime[vers]
	reltimememlimits.update({vers : reltime})

	reltime = float(timetimeouts[vers]) / runtime[vers]
	reltimetimeouts.update({vers : reltime})

	reltime = float(timereaderrors[vers]) / runtime[vers]
	reltimereaderrors.update({vers : reltime})

	reltime = float(timesolved[vers]) / runtime[vers]
	reltimesolved.update({vers : reltime})

faildata = {'fails': reltimefails.values(), 'aborts': reltimeaborts.values(), 'memlimits': reltimememlimits.values(), 
	'timeouts': reltimetimeouts.values(), 'readerrors': reltimereaderrors.values(), 'solved': reltimesolved.values()}
failbars = pd.DataFrame(data=faildata)
failbars.plot(kind='bar', stacked=True, width=0.4)

# label the stacked bars
labelscale = 0.02 
for (ind,row) in failbars.iterrows():
	cumval = 0
	lastleft = True
	for column in failbars:
		val = row.loc[column]
		if not val == 0:
			cumval = cumval + val
			if val < labelscale and not lastleft:
				plt.annotate( round(val,2), xy = (ind+0.2, cumval), horizontalalignment='left', 
					verticalalignment='top', fontsize=6 )
				lastleft = True
			elif val < labelscale:
				plt.annotate( round(val,2), xy = (ind-0.2, cumval), horizontalalignment='right', 
					verticalalignment='top', fontsize=6 )
				lastleft = False
			else:
				plt.annotate( round(val,2), xy = (ind, cumval-0.005), horizontalalignment='center', 
					verticalalignment='top', fontsize=6 )

plt.xticks(range(len(fails)), fails.keys(), rotation=90)
setbarplotparams(1)
plt.ylim(ymax=1.1)
plt.ylabel('Percentage of runtime', size=7)
plt.tick_params(axis='y', which='major', labelsize=7)

ax1 = plt.subplot(111)
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=False, prop={'size': 'small'}, framealpha=1.0)

# if the number of instances differs

if ninstances < 0:
	plt.figtext(.01,.01,'The total number of instances in the test (per version) was unknown or differed.', 
		size='x-small')

plt.savefig(outdir + '/timecomparisonperstatus.pdf')			# name of image
