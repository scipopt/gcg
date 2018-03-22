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

for key in datasets.keys():
	# crop the filenames (keys in datasets) by removing res_ ... .pkl
	croppedkey = key.split('/')[-1].replace('res_', '').replace('.pkl', '')
	fails[croppedkey] = 0
	for status in datasets[key]['status']:
		if status != 'ok':
			fails[croppedkey] = fails[croppedkey] + 1

# Plot results
fig = plt.figure()
plt.title("Number of unsolved instances")
plt.xlabel("GCG Version")
plt.ylim(ymin=0)
plt.bar(range(len(fails)), fails.values(), align='center')
plt.xticks(range(len(fails)), fails.keys())
plt.savefig('images/failcomparison.pdf')
