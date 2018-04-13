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

# Count number of not "ok" instances (failed/aborted/timeout/...) for each res file
fails = {}
highestfails = 0

for key in datasets.keys():
	# crop the filenames (keys in datasets) by removing res_ ... .pkl
	croppedkey = key.split('/')[-1].replace('res_', '').replace('.pkl', '')
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

for item in bars:
        height = item.get_height()
        ax.text(item.get_x()+item.get_width()/2., 1.01*height,
                '%d' % int(height),
                ha='center')

plt.savefig('images/failcomparison.pdf')			# name of image
