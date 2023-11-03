#!/usr/bin/env python3
# Comparison plots for res files.
# Takes res files in pickles format (aquired by parseres.py).
# Python 3

import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import collections
import math
from tikzplotlib import save as tikz_save

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
	outdir = sys.argv[2]
	outdirset = True
	if not os.path.exists(outdir):
	    os.makedirs(outdir)

# Get premade res data
datasets = {}
sumsets = {}
timelimitset = {}

filenames = []
sumnames = []
timelimitnames = []
readmeexists = False
testset = 'short'

comparesettings = False
settingslist = []

orderByCommand = False
sortedbranches = []

# helper function to add an increasing number to duplicate instance names in the ordering
def rename_duplicates( old ):
	seen = {}
	for x in old:
		if x in seen:
			seen[x] += 1
			print(seen[x])
			yield "%s%d" % (x, seen[x])
		else:
			seen[x] = 0
			yield x

for resfile in os.listdir(resdir):
	if resfile.endswith('.pkl') and resfile.startswith('res_'):
		datasets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		filenames.append(resfile)
	elif resfile.endswith('.pkl') and resfile.startswith('sumres_'):
		sumsets[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		sumnames.append(resfile)
	elif resfile.endswith('.pkl') and resfile.startswith('timelimit_'):
		timelimitset[resfile] = pd.read_pickle(os.path.join(resdir, resfile))
		timelimitnames.append(resfile)
	elif resfile.endswith('.txt') and resfile.startswith('readme'):
		# Check for testset name
		filename = resdir + '/' + resfile
		readfile = open(filename, 'r')
		notice = False
		parameterLine = False
		globalLine = False
		globalflags=[]
		for line in readfile:
			if line.startswith('  Testset:'):
				readmeexists = True
				testset = line.split(' ')[-1].split('\n')[0] # line is of form "Testset testsetname"
			if line.startswith('Note'):
				notice = True
			# get the ordering from the readme file
			if orderByCommand and line.startswith('  Branch:'):
				sortedbranches.append(line.split(' ')[-1].split('\n')[0])
				if len(sortedbranches) > len(set(sortedbranches)) and not comparesettings:
					print("You entered the same branch twice. Using settings compare mode.")
					comparesettings = True
			# get settings from readme file
			if line.startswith('  Settings:'):
				settingslist.append(line.split(' ')[-1].split('\n')[0])
		if not notice:
			readfile = open(filename, 'a')
			readfile.write("Note: All plots (apart from \"runtimes\") count the runtime of all fails, aborts, timelimits, memlimits and readerrors as running into the timelimit.")

if comparesettings:
	sortedbranches = list(rename_duplicates(sortedbranches))
print("Using ordering: {}".format(sortedbranches))
if comparesettings:
	print("Using settings: {}".format(settingslist))
#for i in range(len(sortedbranches)):
#	sortedbranches[i] = "res_" + sortedbranches[i] + ".pkl"

# sort names alphabetically
ordereddata_temp = collections.OrderedDict(datasets.items())
orderedsum = collections.OrderedDict(sumsets.items())
orderedtimelimit = collections.OrderedDict(timelimitset.items())

# resort names according to readme file, where the arguments as given in the shell script, are saved
ordereddata = collections.OrderedDict()

# initialize ordereddata randomly if no ordering was given
if len(sortedbranches) == 0:
	sortedbranches = [run.split(".pkl")[0].split("res_")[1] for run in ordereddata_temp.keys()]
	
for k in sortedbranches:
    ordereddata["res_{}.pkl".format(k)] = ordereddata_temp["res_{}.pkl".format(k)]

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
	print('--------------------------------------------------------------------------------------------------')
	print('Warning: Not all tests had the same number of instances.')
	print('Did you enter more than one testset? Did all tested versions have access to all testset instances?')
	print('--------------------------------------------------------------------------------------------------')
	ninstances = -1

# Sanity check: check for fails, let the dev know
for res in sumnames:
	if not sumsets[res]['Fail'] == '0':
		print('--------------------------------------------------------------------------------------------------')
		print('Warning: There were some failed runs in the tests. This might influence the significance of the')
		print('comparisons! Recommendation: Check for memlimits, aborts, fails etc. in the tested GCG versions.')
		print('--------------------------------------------------------------------------------------------------')
		break

# Sanity check: check whether the timelimits were indentical for all versions
defaulttimelimit = -1
printwarning = False
for res in timelimitnames:
	if int(defaulttimelimit) == -1:
		defaulttimelimit = orderedtimelimit[res]['timelimit']
	else:
		if int(defaulttimelimit) != int(orderedtimelimit[res]['timelimit']):
			printwarning = True
if printwarning == True:
	print('--------------------------------------------------------------------------------------------------')
	print('Warning: The timelimit of the versions differed. Some plots use the timelimit as a default for all')
	print('fails. Recommendation: Rurun the tests with a global timelimit.')
	print('--------------------------------------------------------------------------------------------------')

# -------------------------------------------------------------------------------------------------------------------------
# Add function to crop filenames
# -------------------------------------------------------------------------------------------------------------------------

def cropkeypkl(key, keyprefix, addlinebreak=False):
	# crop the filenames by removing prefix ... .pkl and add linebreak for very long keys
	croppedkey = key.split('/')[-1].replace(keyprefix, '').replace('.pkl', '')
	if len(croppedkey) > maxstringlen:
		charlist = list(croppedkey)
		ninserts = len(charlist)/maxstringlen
		while ninserts > 0:
			if addlinebreak:
				charlist.insert(int(ninserts*maxstringlen), '\n')
			ninserts = ninserts - 1
		croppedkey = ''.join(charlist)
	return croppedkey

# -------------------------------------------------------------------------------------------------------------------------
# Get some statistics for each res file to be used in the plots
# -------------------------------------------------------------------------------------------------------------------------

maxstringlen = 12

versions = []
timelimits = {}

fails = {}
aborts = {}
memlimits = {}
timeouts = {}

timefails = {}
timeaborts = {}
timememlimits = {}
timetimeouts = {}
timesolved = {}

nsolved = {}

timeperinstance = {}

highestfails = 0
tempruntime = {}
highesttime = 0

# extract timelimits
for key in list(orderedtimelimit.keys()):
	croppedkey = cropkeypkl(key, 'timelimit_')
	timelimits[croppedkey] = int(orderedtimelimit[key]['timelimit'])

for key in list(ordereddata.keys()):
	# crop the filenames (keys in ordereddata) by removing res_ ... .pkl and add linebreak for very long keys
	croppedkey = cropkeypkl(key, 'res_')
	versions.append(croppedkey)

	# get fail types and their amounts
	fails[croppedkey] = 0
	aborts[croppedkey] = 0
	memlimits[croppedkey] = 0
	timeouts[croppedkey] = 0
	for status in ordereddata[key]['status']:
		if status == 'fail' or status == 'readerror':
			fails[croppedkey] = fails[croppedkey] + 1
		elif status == 'abort':
			aborts[croppedkey] = aborts[croppedkey] + 1
		elif status == 'memlimit':
			memlimits[croppedkey] = memlimits[croppedkey] + 1
		elif status == 'timeout':
			timeouts[croppedkey] = timeouts[croppedkey] + 1

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

	# get runtime per instance for each version
	temptimeperinstance = {}
	for i in range(ninstances):
		tempinsname = ordereddata[key]['Name'][i]
		# the instance names might not be unique but they will appear in the same order in all versions
		if tempinsname in temptimeperinstance:
			while tempinsname in temptimeperinstance:
				tempinsname = tempinsname + '_'
		if ordereddata[key]['status'][i] == 'fail' or ordereddata[key]['status'][i] == 'readerror' or ordereddata[key]['status'][i] == 'abort' or ordereddata[key]['status'][i] == 'abort'or ordereddata[key]['status'][i] == 'timeout':
			temptimeperinstance.update({tempinsname: timelimits[croppedkey]})
		else:
			temptimeperinstance.update({tempinsname: ordereddata[key]['TotalTime'][i]})

	timeperinstance.update({croppedkey: temptimeperinstance})

	# get runtime per status
	timefails[croppedkey] = 0
	timeaborts[croppedkey] = 0
	timememlimits[croppedkey] = 0
	timetimeouts[croppedkey] = 0
	timesolved[croppedkey] = 0
	nsolved[croppedkey] = 0
	tablelength = len(ordereddata[key].index)
	for i in range(0, tablelength-1):
		if ordereddata[key]['status'][i] == 'fail':
			timefails[croppedkey] = timefails[croppedkey] + timelimits[croppedkey]
		elif ordereddata[key]['status'][i] == 'readerror':
			timefails[croppedkey] = timefails[croppedkey] + timelimits[croppedkey]
		elif ordereddata[key]['status'][i] == 'abort':
			timeaborts[croppedkey] = timeaborts[croppedkey] + timelimits[croppedkey]
		elif ordereddata[key]['status'][i] == 'abort':
			timememlimits[croppedkey] = timememlimits[croppedkey] + timelimits[croppedkey]
		elif ordereddata[key]['status'][i] == 'timeout':
			timetimeouts[croppedkey] = timetimeouts[croppedkey] + float(ordereddata[key]['TotalTime'][i])
		else:
			timesolved[croppedkey] = timesolved[croppedkey] + float(ordereddata[key]['TotalTime'][i])
			nsolved[croppedkey] = nsolved[croppedkey] + 1

	# round up runtime per status
	timefails[croppedkey] = math.ceil(timefails[croppedkey])
	timeaborts[croppedkey] = math.ceil(timeaborts[croppedkey])
	timememlimits[croppedkey] = math.ceil(timememlimits[croppedkey])
	timetimeouts[croppedkey] = math.ceil(timetimeouts[croppedkey])
	timesolved[croppedkey] = math.ceil(timesolved[croppedkey])

# DO NOT order statistics by keys
nversions = len(versions)
versions = versions

fails = collections.OrderedDict(sorted(fails.items()))
aborts = collections.OrderedDict(sorted(aborts.items()))
memlimits = collections.OrderedDict(sorted(memlimits.items()))
timeouts = collections.OrderedDict(sorted(timeouts.items()))
runtime_temp = collections.OrderedDict(sorted(tempruntime.items()))###

# resort names according to readme file, where the arguments as given in the shell script, are saved
runtime = collections.OrderedDict()
for k in sortedbranches:
	runtime[k] = runtime_temp[k]


timefails = collections.OrderedDict(sorted(timefails.items()))
timeaborts = collections.OrderedDict(sorted(timeaborts.items()))
timememlimits = collections.OrderedDict(sorted(timememlimits.items()))
timetimeouts = collections.OrderedDict(sorted(timetimeouts.items()))
timesolved = collections.OrderedDict(sorted(timesolved.items()))

# add a runtime where every fail type is counted as timelimit
totalruntime = collections.OrderedDict()
highesttotalruntime = 0.0
for key in list(runtime.keys()):
	time = 0.0
	for instance in list(timeperinstance[key].keys()):
		time = time + float(timeperinstance[key][instance])
	totalruntime.update({key: time})
	highesttotalruntime = max(highesttotalruntime, time)

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
	return

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
	plt.tick_params(axis='x', labelsize=7)
	return

# -------------------------------------------------------------------------------------------------------------------------
# 1) Plot how many instances were unsolved per version
# -------------------------------------------------------------------------------------------------------------------------

fig = plt.figure()
ax = plt.axes()
plt.title('Number of unsolved instances')
if not comparesettings:
	plt.xlabel('GCG Version')
else:
	plt.xlabel('GCG Settings')

faildata = collections.OrderedDict([('aborts', list(aborts.values())), ('fails', list(fails.values())), ('memlimits', list(memlimits.values())),
	('timeouts', list(timeouts.values()))])
failbars = pd.DataFrame(data=faildata)
failbars.plot(kind='bar', stacked=True)

# calculate highest bar length
barheight = 0
for vers in list(timefails.keys()):
	barheight = max(barheight, aborts[vers] + fails[vers] + memlimits[vers] + timeouts[vers])

barheight = barheight + .1*barheight

# label the stacked bars
labelscale = 0.02*barheight
for (ind,row) in failbars.iterrows():
	cumval = 0
	lastleft = True
	for column in failbars:
		val = row.loc[column]
		if not val == 0:
			cumval = cumval + val
			if val < labelscale and not lastleft:
				plt.annotate( val, xy = (ind-.3, cumval - .5), horizontalalignment='left', verticalalignment='top',
					fontsize=6 )
				lastleft = True
			elif val < labelscale:
				plt.annotate( val, xy = (ind+.3, cumval - .5), horizontalalignment='right', verticalalignment='top',
					fontsize=6 )
				lastleft = False
			else:
				plt.annotate( val, xy = (ind, cumval - .5), horizontalalignment='center', verticalalignment='top',
					fontsize=6 )

if not comparesettings:
	plt.xticks(list(range(len(fails))), list([cropkeypkl(key,"",True) for key in fails.keys()]), rotation=90)
else:
	plt.xticks(list(range(len(fails))), list(settingslist), rotation=90)
setbarplotparams(int(float(highestfails)))

ax1 = plt.subplot(1,1,1,label="1")
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5, fancybox=False, prop={'size': 'small'}, framealpha=1.0)

# if the number of instances differs
stringninstances = 'unknown or differed'
if ninstances >= 0:
	stringninstances = str(ninstances)

plt.figtext(.01,.01,'The total number of instances in the test (per version) was ' + stringninstances, size='x-small')
if comparesettings:
	plt.figtext(.01,.060,'Branch: ' + list(sortedbranches)[0], size='x-small')
plt.figtext(.01,.035,'Testset: ' + testset, size='x-small')
plt.subplots_adjust(bottom=0.2)

plt.savefig(outdir + '/failcomparison.pdf')			# name of image
tikz_save(outdir + '/failcomparison.tikz',
	axis_height = '\\figureheight',
	axis_width= '\\figurewidth')

# -------------------------------------------------------------------------------------------------------------------------
# 2) Plot runtime per version
# -------------------------------------------------------------------------------------------------------------------------

fig = plt.figure()
ax = plt.axes()
plt.title('Runtime per version')
plt.ylabel('Runtime in seconds')

bars = plt.bar(list(range(len(runtime))), list(runtime.values()), align='center')
if not comparesettings:
	plt.xticks(list(range(len(runtime))), list([cropkeypkl(key,"",True) for key in runtime.keys()]), rotation=90)
else:
	plt.xticks(list(range(len(runtime))), list(settingslist), rotation=90)
setbarplotparams(highesttime)
labelbars(bars, highesttime)

if comparesettings:
	plt.figtext(.01,.035,'Branch: ' + list(sortedbranches)[0], size='x-small')
plt.figtext(.01,.01,'Testset: ' + testset, size='x-small')

plt.savefig(outdir + '/runtimes.pdf')				# name of image
tikz_save(outdir + '/runtimes.tikz',
	axis_height = '\\figureheight',
	axis_width= '\\figurewidth')

# -------------------------------------------------------------------------------------------------------------------------
# 3.a) Some helper functions for cumulative time differences
# -------------------------------------------------------------------------------------------------------------------------

# sums up runtimes of different instances given a list of names and a dict with (name, time) tuples
# (avoid naming errors by checking whether name exists)
def sumruntimes(namelist, instimelist):
	res = 0.0
	for insname in namelist:
		if insname in instimelist:
			res = res + float(instimelist[insname])
	return res

# calculate speedup factor given a list of (version, value) tuples and an index
def calcspeedup(vallist, i):
	assert i > 0
	if not vallist[i][1] == 0:
		speedup = float(vallist[i-1][1]) / vallist[i][1]
	else:
		speedup = 0
	return speedup

# add value to former one: get a (key, value) dict, a key with index > 0 and the current difference
def addtoformer(valuedict, key, diff):
	cumitems = list(valuedict.items())
	assert len(cumitems)-1 >= 0
	cumsum = cumitems[len(cumitems)-1][1]

	res = float(diff) * cumsum
	return res

# -------------------------------------------------------------------------------------------------------------------------
# 3) Plot runtime comparison
# -------------------------------------------------------------------------------------------------------------------------

# Calculate version-to-version speedup

items = list(totalruntime.items())
assert(len(items) == nversions)
if nversions < 2:
	print('Enter more than one GCG version to generate a runtime comparison plot.')
else:
	# get instance names that originally (in first version) ran in under 10, 100, 1000 seconds
	names10 = []
	names100 = []
	names1000 = []
	nameslong = []

	# only compute intervals if > 2 versions
	if nversions > 2:
		(firstvers, instances) = list(timeperinstance.items())[len(list(timeperinstance.items()))-1]
		# get names of instances running in certain intervals on the latest version
		for instancename in list(instances.keys()):
			if float(instances[instancename]) < 10.0:
				names10.append(instancename)
			elif float(instances[instancename]) < 100.0:
				names100.append(instancename)
			elif float(instances[instancename]) < 1000.0:
				names1000.append(instancename)
			else:
				nameslong.append(instancename)

	# get sum of runtimes of these instances
	runtimes10 = collections.OrderedDict()
	runtimes100 = collections.OrderedDict()
	runtimes1000 = collections.OrderedDict()
	runtimeslong = collections.OrderedDict()

	if nversions > 2:
		for vers in list(timeperinstance.keys()):
			runtimes10[vers] = sumruntimes(names10, timeperinstance[vers])
			runtimes100[vers] = sumruntimes(names100, timeperinstance[vers])
			runtimes1000[vers] = sumruntimes(names1000, timeperinstance[vers])
			runtimeslong[vers] = sumruntimes(nameslong, timeperinstance[vers])

		# convert the runtimes for easier access
		runtimes10 = sorted(list(runtimes10.items()))
		runtimes100 = sorted(list(runtimes100.items()))
		runtimes1000 = sorted(list(runtimes1000.items()))
		runtimeslong = sorted(list(runtimeslong.items()))

	# prepare variables
	highestdiff = 0
	lowestdiff = 0

	runtimecomp = collections.OrderedDict()
	cumulative = collections.OrderedDict() # overall cumulative speedup
	cum10 = collections.OrderedDict() # cumulative speedup for instances with original runtime <10s
	cum100 = collections.OrderedDict() # cumulative speedup for instances with original runtime <100s
	cum1000 = collections.OrderedDict() # cumulative speedup for instances with original runtime <1000s
	cumlong = collections.OrderedDict() # cumulative speedup for instances with original runtime >1000s

	highestcum = 0
	lowestcum = 0
	highestcum10 = 0
	lowestcum10 = 0
	highestcum100 = 0
	lowestcum100 = 0
	highestcum1000 = 0
	lowestcum1000 = 0
	highestcumlong = 0
	lowestcumlong = 0

	# calculate the different times per version
	for i in range(nversions):
		diff = 0
		diff10 = 0
		diff100 = 0
		diff1000 = 0
		difflong = 0

		if i > 0:
			# from the second item on calculate the version speed differences (speedup)
			name = items[i-1][0] + '\n->\n' + items[i][0]
			if comparesettings:
				name = settingslist[i-1] + ' -> ' + settingslist[i]
			diff = calcspeedup(items, i)
			runtimecomp[name] = diff
			highestdiff = max(float(diff), float(highestdiff))
			lowestdiff = min(float(diff), float(lowestdiff))

			if nversions > 2:
				diff10 = calcspeedup(runtimes10, i)
				diff100 = calcspeedup(runtimes100, i)
				diff1000 = calcspeedup(runtimes1000, i)
				difflong = calcspeedup(runtimeslong, i)

			# for the first one set initial cumulative difference values
			if i == 1:
				cumulative[name] = highestcum = lowestcum = diff
				if nversions > 2:
					cum10[name] = highestcum10 = lowestcum10 = diff10
					cum100[name] = highestcum100 = lowestcum100 = diff100
					cum1000[name] = highestcum1000 = lowestcum1000 = diff1000
					cumlong[name] = highestcumcumlong = lowestcumcumlong = difflong

			# for all following add the last value to current diff
			else:
				cumulative[name] = addtoformer(cumulative, name, diff)
				highestcum = max(float(highestcum), float(cumulative[name]))
				lowestcum = min(float(lowestcum), float(cumulative[name]))

				if nversions > 2:
					cum10[name] = addtoformer(cum10, name, diff10)
					highestcum10 = max(float(highestcum10), float(cum10[name]))
					lowestcum10 = min(float(lowestcum10), float(cum10[name]))

					cum100[name] = addtoformer(cum100, name, diff100)
					highestcum100 = max(float(highestcum100), float(cum100[name]))
					lowestcum100 = min(float(lowestcum100), float(cum100[name]))

					cum1000[name] = addtoformer(cum1000, name, diff1000)
					highestcum1000 = max(float(highestcum1000), float(cum1000[name]))
					lowestcum1000 = min(float(lowestcum1000), float(cum1000[name]))

					cumlong[name] = addtoformer(cumlong, name, difflong)
					highestcumlong = max(float(highestcumlong), float(cumlong[name]))
					lowestcumlong = min(float(lowestcumlong), float(cumlong[name]))

	#determine axis min/max
	if nversions > 2:
		axmin = min(lowestcum, lowestdiff, lowestcum10, lowestcum100, lowestcum1000, lowestcumlong)
		axmax = max(highestcum, highestdiff, highestcum10, highestcum100, highestcum1000, highestcumlong)
	else:
		axmin = min(lowestcum, lowestdiff)
		axmax = max(highestcum, highestdiff)

	# make space for bar labels
	longestbar = max(axmax, abs(axmin))
	axmin = axmin - 0.1*longestbar
	axmax = axmax + 0.1*longestbar

	# first plot version-to-version comparison bars
	fig, ax1 = plt.subplots()
	bar1 = ax1.bar(list(range(len(runtimecomp))), list(runtimecomp.values()), color='b')
	if not comparesettings:
		plt.xticks(list(range(len(runtimecomp))), list([cropkeypkl(key,"",True) for key in runtimecomp.keys()]), rotation=90)
	else:
		plt.xticks(list(range(len(runtimecomp))), list(settingslist), rotation=90)
	plt.tick_params(axis='x', labelsize=5)
	ax1.set_ylabel('Speedup factor', color='b')
	ax1.tick_params('y', colors='b')
	ax1.set_ylim(ymin=axmin, ymax=axmax)

	#longestbar = max(highestdiff, abs(lowestdiff))
	labelbars(bar1, longestbar)

	# plot cumulative speedup if there is more than one bar
	nbars = list(range(len(runtimecomp)))
	if len(items) > 2:
		ax2 = ax1.twinx()
		ax2.set_ylabel('Cumulative speedup factor', color='r')
		ax2.tick_params('y', colors='r')
		ax2.set_ylim(ymin=axmin, ymax=axmax)

		ax2.plot(nbars, list(cumulative.values()), 'r-', label='overall')
		ax2.axhline(y=0, color='xkcd:slate')

		if nversions > 2:
			ax2.plot(nbars, list(cum10.values()), 'xkcd:light orange', label='<10s')
			ax2.plot(nbars, list(cum100.values()), 'xkcd:orange', label='[10,100)s')
			ax2.plot(nbars, list(cum1000.values()), 'xkcd:dark orange', label='[100,1000)s')
			ax2.plot(nbars, list(cumlong.values()), 'xkcd:reddy brown', label='>1000s')

		ax2.legend(loc='upper right', prop={'size': 'x-small'})

	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	plt.title('Version-to-version runtime comparison')

	stringninstances = 'unknown or differed'
	if ninstances >= 0:
		stringninstances = str(ninstances)

	if nversions > 2:
		plt.figtext(.01,0,'The total number of instances in the test (per version) was ' + stringninstances + '.\n' +
			'Amount of instances running in the latest version in: \n<10s: ' + str(len(names10)) + '.\n' +
			'[10,100)s: ' + str(len(names100)) + '.\n' +
			'[100,1000)s: ' + str(len(names1000)) + '.\n' +
			'>1000s: ' + str(len(nameslong)), size='x-small')
		if comparesettings:
			plt.figtext(.06,.06,'Settings: ' + list(runtimecomp.keys())[0], size='x-small')
			plt.figtext(.01,.06,'Branch: ' + list(sortedbranches)[0], size='x-small')
		plt.figtext(.01,.035,'Testset: ' + testset, size='x-small')
		plt.subplots_adjust(bottom=0.25)
	else:
		plt.figtext(.01,.01,'The total number of instances in the test (per version) was ' + stringninstances + '.'
			,size='x-small')
		if comparesettings:
			plt.figtext(.60,.035,'Settings: ' + list(runtimecomp.keys())[0], size='x-small')
			plt.figtext(.01,.06,'Branch: ' + list(sortedbranches)[0], size='x-small')
		plt.figtext(.01,.035,'Testset: ' + testset, size='x-small')

	plt.savefig(outdir + '/runtimecomparison.pdf')			# name of image
	tikz_save(outdir + '/runtimecomparison.tikz',
		axis_height = '\\figureheight',
		axis_width= '\\figurewidth')

# -------------------------------------------------------------------------------------------------------------------------
# 4) Plot time per status category (fail categories and solved)
# -------------------------------------------------------------------------------------------------------------------------

fig = plt.figure()
ax = plt.axes()

faildata = collections.OrderedDict([('aborts', list(timeaborts.values())), ('fails', list(timefails.values())), ('memlimits', list(timememlimits.values())),
	('timeouts', list(timetimeouts.values())), ('solved',list(timesolved.values()))])
failbars = pd.DataFrame(data=faildata)
failbars.plot(kind='bar', stacked=True, width=0.4)

# calculate highest bar length
barheight = highesttime + .1*highesttime
for vers in list(timefails.keys()):
	barheight = max(barheight, timefails[vers] + timeaborts[vers] + timememlimits[vers] + timetimeouts[vers] + timesolved[vers])

barheight = barheight + .1*barheight

# label the stacked bars
labelscale = 0.02*barheight
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

if not comparesettings:
	plt.xticks(list(range(len(fails))), list([cropkeypkl(key,"",True) for key in fails.keys()]), rotation=90)
else:
	plt.xticks(list(range(len(fails))), list(settingslist), rotation=90)
setbarplotparams(1)
plt.ylim(ymax=barheight)
plt.ylabel('Runtime in seconds', size=7)
plt.tick_params(axis='y', labelsize=7)

ax1 = plt.subplot(1,1,1,label="2")
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=False, prop={'size': 'small'}, framealpha=1.0)

# if the number of instances differs

if ninstances < 0:
	plt.figtext(.01,.01,'The total number of instances in the test (per version) was unknown or differed.' +
	'Testset: ' + testset, size='x-small')
else:
	if comparesettings:
		plt.figtext(.01,.035,'Branch: ' + list(sortedbranches)[0], size='x-small')
	plt.figtext(.01,.01,'Testset: ' + testset, size='x-small')

plt.subplots_adjust(left=0.1)

plt.savefig(outdir + '/timecomparisonperstatus.pdf')			# name of image
tikz_save(outdir + '/timecomparisonperstatus.tikz',
           axis_height = '\\figureheight',
           axis_width= '\\figurewidth')

# -------------------------------------------------------------------------------------------------------------------------
# 5) Plot average runtime of solved instances per version
# -------------------------------------------------------------------------------------------------------------------------

fig = plt.figure()
ax = plt.axes()
plt.title('Average runtime of solved instances')
if not comparesettings:
	plt.xlabel('GCG Version')
else:
	plt.xlabel('GCG Settings')
plt.ylabel('Average runtime in seconds')

avsolved = collections.OrderedDict()
highestavsolved = 0

for vers in versions:
	if not nsolved[vers] == 0:
		avtime = float(timesolved[vers]) / nsolved[vers]
	else:
		avtime = timelimits[vers]
	avsolved.update({vers : avtime})
	highestavsolved = max(highestavsolved, avtime)

bars = plt.bar(list(range(len(avsolved))), list(avsolved.values()), align='center')
if not comparesettings:
	plt.xticks(list(range(len(avsolved))), list([cropkeypkl(key,"",True) for key in avsolved.keys()]), rotation=90)
else:
	plt.xticks(list(range(len(avsolved))), list(settingslist), rotation=90)
setbarplotparams(highestavsolved)
labelbars(bars, highestavsolved)
if comparesettings:
	plt.figtext(.01,.035,'Branch: ' + list(sortedbranches)[0], size='x-small')
plt.figtext(.01,.01,'Testset: ' + testset, size='x-small')

plt.savefig(outdir + '/averagesolvetime.pdf')				# name of image
tikz_save(outdir + '/averagesolvetime.tikz',
           axis_height = '\\figureheight',
           axis_width= '\\figurewidth')
