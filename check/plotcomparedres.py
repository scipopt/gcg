#!/usr/bin/env python
# Comparison plots for res files. 
# Takes res files in pickles format (aquired by parseres.py).
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# 1) Plot how many instances were unsolved per version

# Get premade res data
datasets = {}
filenames = []

for file in os.listdir("pickles/"):
	if file.endswith(".pkl") and file.startswith("res_"):
		datasets[file] = pd.read_pickle(os.path.join("pickles/", file))
		filenames.append(file)

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
	# countfails
	fails[croppedkey] = 0
	for status in datasets[key]['status']:
		if status != 'ok' and status != 'solved' and status != 'solvednotverified':
			fails[croppedkey] = fails[croppedkey] + 1
			if highestfails < fails[croppedkey]:
				highestfails = fails[croppedkey]

# Plot results
fig = plt.figure()
ax = plt.axes()        
plt.title("Number of unsolved instances")
plt.xlabel("GCG Version")
plt.ylim(ymin=0)
plt.ylim(ymax=highestfails+(highestfails/10))			# max y value is set to more than highest value
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

plt.savefig('images/failcomparison.pdf')			# name of image
