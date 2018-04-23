#!/usr/bin/env python
# This script reads *.out files from whole testsets created by GCG's make test
# and parses its statistics into a pandas dataframe. 
# Very hacky, beware!
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# check command line arguments
if len(sys.argv) < 2:
	sys.exit("Usage: ./parseout.py OUTFILE")

# help functions
def ct(string):
	sim = string.strip()[1:-1]
	return sim.upper()

def real_path(string):
	splitted = string.split('/')
	check = False
	ret = ""
	for w in splitted:
		if w == "check":
			check = True
		elif check:
			ret += "/" + w
		else:
			continue
	return ret[1:]

# data dictionary
data = {
	'OVERALL TIME': [],
	'STATUS': [],
	'ROOT NODE TIME': [],
	'DUAL BOUNDS': [],
	'HEURISTICS TIME': [],
	'HEURISTICS CALLS': [],
	'HEURISTICS FOUND': [],
	'CUTS TIME': [],
	'CUTS CALLS': [],
	'CUTS FOUND': [],
	'CUTS APPLIED': [],
	'FARKAS TIME': [],
	'MASTER TIME': [],
	'PRICING TIME': [],
	'MIP OR KNAPSACK': [],
	'DEGENERACY': [],
	'CONS LINEAR': [],
	'CONS KNAPSACK': [],
	'CONS LOGICOR': [],
	'CONS SETPPC': [],
	'CONS VARBOUND': [],
	'LINKING VARS': [],
	'NBLOCKS': [],
	'AGGREGATIONS': [],
	'SOLUTIONS FOUND': [],
	'SOLUTIONS FIRST': [],
	'SOLUTIONS BEST': [],
	'PD INTEGRAL': [],
	'MASTER NCONSS': [],
	'MASTER NVARS': [],
	'BNB TREE NODES': [],
	'BNB TREE LEFT': [],
	'BNB TREE DEPTH': [],
	'LP CALLS': [],
	'LP TIME': [],
	'LP ITERATIONS': [],
	'LP FILE': [],
	'DEC FILE': [],
}


# instance names
index = []

# checkout outfile
outfile = sys.argv[1]
fh = open(outfile, 'r')

# initialize variables
search = ""
it = 0
opstat = False
ot = False
status = False

# write in dictionary
for line in fh:
	# get instance name
	if line.startswith("read problem") and '.dec' not in line and '.blk' not in line:
		index.append(line.split()[2].split('/')[-1].split('.')[0])

	# get instance lp
	if line.startswith("read problem") and '.dec' not in line and '.blk' not in line:
		data['LP FILE'].append(real_path(line.split()[2][1:-1]))
	
	# get instance dec
	if line.startswith("read problem") and ('.dec' in line or '.blk' in line):
		data['DEC FILE'].append(real_path(line.split()[2][1:-1]))
	
	# reading of master stats finished
	if line.startswith("Original Program statistics:"):
		opstat = True
		continue

	# get overall time
	if line.startswith("Solving Time (sec)") and not ot:
		data['OVERALL TIME'].append(float(line.split(':')[1]))
		ot = True
		continue

	# get status
	if line.startswith("SCIP Status") and not status:
		if line.split(':')[1].strip() == "problem is solved [optimal solution found]":
			data['STATUS'].append(1)
		elif line.split(':')[1].strip() == "problem is solved [infeasible]":
			data['STATUS'].append(2)
		elif line.split(':')[1].strip() == "solving was interrupted [time limit reached]":
			data['STATUS'].append(3)
		elif line.split(':')[1].strip() == "solving was interrupted [memory limit reached]":
			data['STATUS'].append(4)
		elif line.split(':')[1].strip() == "solving was interrupted [node limit reached]":
			data['STATUS'].append(5)
		else:
			data['STATUS'].append(0)
		status = True

	# get root node time
	if line.startswith("Time in root node:"):
		data['ROOT NODE TIME'].append(float(line.split(':')[1]))
		continue

	# get degeneracy
	if line.startswith("Degeneracy:"):
		search = "DEGENERACY"
		data['DEGENERACY'].append([])
		continue

	if search == "DEGENERACY":
		if line.startswith("Dual Bounds:"):
			search = "DUALS"	# no empty string here!
			data['DUAL BOUNDS'].append([])
			continue
		data['DEGENERACY'][-1].append((int(line.split(':')[0]), float(line.split(':')[1])))
		continue

	# get dual bound development
	if line.startswith("Dual Bounds:"):
		search = "DUALS"
		data['DUAL BOUNDS'].append([])
		continue

	if search == "DUALS":
		if line.startswith("GCG"):
			search = ""
			continue
		data['DUAL BOUNDS'][-1].append((int(line.split(':')[0]), float(line.split(':')[1])))
		continue

	# get successful heuristics
	if line.startswith("Primal Heuristics"):
		search = "HEURISTICS"
		data['HEURISTICS TIME'].append(0.)
		data['HEURISTICS CALLS'].append(0)
		data['HEURISTICS FOUND'].append(0)
		continue

	if search == "HEURISTICS":
		if not line.startswith("Diving Statistics"):
			if line.split(':')[1].split()[0].replace('.', '', 1).isdigit():
				data['HEURISTICS TIME'][-1] += float(line.split(':')[1].split()[0])
				if line.split(':')[1].split()[1].replace('.', '', 1).isdigit():
					data['HEURISTICS TIME'][-1] += float(line.split(':')[1].split()[1])
				if line.split(':')[1].split()[2].isdigit():
					data['HEURISTICS CALLS'][-1] += int(line.split(':')[1].split()[2])
				if line.split(':')[1].split()[3].isdigit():
					data['HEURISTICS FOUND'][-1] += int(line.split(':')[1].split()[3])
			else:
				search = ""

	# get cutting plane statistics
	if line.startswith("Separators"):
		search = "CUTS"
		data['CUTS TIME'].append(0.)
		data['CUTS CALLS'].append(0)
		data['CUTS FOUND'].append(0)
		data['CUTS APPLIED'].append(0)
		continue

	if search == "CUTS":
		if not line.startswith("Pricers"):
			offset = 0 if line.split(':')[0].strip() != "cut pool" else -1
			data['CUTS TIME'][-1] += float(line.split(':')[1].split()[0])
			data['CUTS CALLS'][-1] += int(line.split(':')[1].split()[2+offset])
			data['CUTS FOUND'][-1] += int(line.split(':')[1].split()[5+offset])
			if line.split(':')[1].split()[6+offset].isdigit():
				data['CUTS APPLIED'][-1] += int(line.split(':')[1].split()[6+offset])
		else:
			search = ""

	# get Farkas Time and type of Pricing (MIP / Knapsack)
	if line.startswith("Pricing Solver"):
		search = "PRICING"
		data['MIP OR KNAPSACK'].append(0)
		data['FARKAS TIME'].append(0.)
		data['PRICING TIME'].append(0.)
		continue

	if search == "PRICING":
		if line.lstrip().startswith("knapsack"):
			if sum([int(x) for x in line.split(':')[1].split()[:4]]) > 0:
				data['MIP OR KNAPSACK'][-1] += 1
				data['FARKAS TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:6]])
				data['PRICING TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:]])
		elif line.lstrip().startswith("mip"):
			if sum([int(x) for x in line.split(':')[1].split()[:4]]) > 0:
				data['MIP OR KNAPSACK'][-1] += 2
				data['FARKAS TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:6]])
				data['PRICING TIME'][-1] += sum([float(x) for x in line.split(':')[1].split()[4:]])
		else:
			search = ""

	# get Master time
	if line.startswith("Master Program statistics:"):
		search = "MASTER"

	if search == "MASTER" and line.split(':')[0].strip() == "solving":
		data['MASTER TIME'].append(float(line.split(':')[1]))
		search = ""

	# get constraints
	if line.startswith("presolved problem has"):
		search = "CONSS"
		data['CONS LINEAR'].append(0)
		data['CONS KNAPSACK'].append(0)
		data['CONS LOGICOR'].append(0)
		data['CONS SETPPC'].append(0)
		data['CONS VARBOUND'].append(0)
		continue

	if search == "CONSS":
		res = re.search("constraints of type", line)
		if res:
			constype = ct(line[res.end():-1])
			data['CONS ' + constype][-1] = int(line[:res.start()])
		else:
			search = ""

	# get number of blocks
	if line.startswith("Decomp statistics"):
		search = "BLOCKS"

	if search == "BLOCKS":
		if line.lstrip().startswith("blocks"):
			data['NBLOCKS'].append(int(line.split(':')[1]))
		if line.lstrip().startswith("aggr. blocks"):
			data['AGGREGATIONS'].append(int(line.split(':')[1]))
			search = ""

	# get solution statistics
	if line.startswith("Solution") and not opstat:
		search = "SOLUTION"
		continue

	if search == "SOLUTION":
		if line.lstrip().startswith("Solutions found"):
			data['SOLUTIONS FOUND'].append(int(line.split(':')[1].split()[0]))
		elif line.lstrip().startswith("First Solution"):
			data['SOLUTIONS FIRST'].append(float(line.split(':')[1].split()[7]))
		elif line.lstrip().startswith("Primal Bound") and data['SOLUTIONS FOUND'][-1] > 0:
			data['SOLUTIONS BEST'].append(float(line.split(':')[1].split()[7]))
		elif line.lstrip().startswith("Avg. Gap") and data['SOLUTIONS FOUND'][-1] > 0:
			data['PD INTEGRAL'].append(float(line.split(':')[1].split('%')[1].split()[0][1:]))
			search = ""
		else:
			continue
	
	# get master statistics
	if line.startswith("Master statistics"):
		search = "MASTER STATS"
		continue

	if search == "MASTER STATS":
		if line.lstrip().startswith("master"):
			data['MASTER NCONSS'].append(int(line.split(':')[1].split()[5]))
			data['MASTER NVARS'].append(int(line.split(':')[1].split()[0]))

	if line.startswith("Number of LinkingVars:"):
		data['LINKING VARS'].append(int(line.split(':')[1]))

	# get Branch-and-Bound Tree stats
	if line.startswith("B&B Tree") and opstat:
		search = "BNB"
		continue

	if search == "BNB":
		if line.lstrip().startswith("nodes (total)"):
			data['BNB TREE NODES'].append(int(line.split(':')[1].split()[0]))
		elif line.lstrip().startswith("nodes left"):
			data['BNB TREE LEFT'].append(int(line.split(':')[1]))
		elif line.lstrip().startswith("max depth (total)"):
			data['BNB TREE DEPTH'].append(int(line.split(':')[1]))
	
	# get LP stats
	if line.startswith("LP") and not opstat:
		search = "LP"
		continue

	if search == "LP":
		if line.lstrip().startswith("primal LP"):
			data['LP CALLS'].append(int(line.split(':')[1].split()[1]))
			data['LP TIME'].append(float(line.split(':')[1].split()[0]))
			data['LP ITERATIONS'].append(int(line.split(':')[1].split()[2]))
		elif line.lstrip().startswith("dual LP"):
			data['LP CALLS'][-1] += int(line.split(':')[1].split()[1])
			data['LP TIME'][-1] += float(line.split(':')[1].split()[0])
			data['LP ITERATIONS'][-1] += int(line.split(':')[1].split()[2])


	# sync point
	if line.startswith("=ready="):
		it += 1
		search = ""
		opstat = False
		ot = False
		status = False
		if len(data['OVERALL TIME']) < it:
			data['OVERALL TIME'].append(float('NaN')) 
		if len(data['STATUS']) < it:
			data['STATUS'].append(-1)
		if len(data['ROOT NODE TIME']) < it:
			data['ROOT NODE TIME'].append(float('NaN')) 

		if len(data['DUAL BOUNDS']) < it:
			data['DUAL BOUNDS'].append([float('NaN')]) 

		if len(data['HEURISTICS TIME']) < it:
			data['HEURISTICS TIME'].append(float('NaN'))
		elif len(data['HEURISTICS TIME']) == it + 1:
			data['HEURISTICS TIME'][-2] += data['HEURISTICS TIME'][-1]
			data['HEURISTICS TIME'] = data['HEURISTICS TIME'][:-1]
		if len(data['HEURISTICS CALLS']) < it:
			data['HEURISTICS CALLS'].append(-1)
		elif len(data['HEURISTICS CALLS']) == it + 1:
			data['HEURISTICS CALLS'][-2] += data['HEURISTICS CALLS'][-1]
			data['HEURISTICS CALLS'] = data['HEURISTICS CALLS'][:-1]
		if len(data['HEURISTICS FOUND']) < it:
			data['HEURISTICS FOUND'].append(-1)
		elif len(data['HEURISTICS FOUND']) == it + 1:
			data['HEURISTICS FOUND'][-2] += data['HEURISTICS FOUND'][-1]
			data['HEURISTICS FOUND'] = data['HEURISTICS FOUND'][:-1]

		if len(data['CUTS TIME']) < it:
			data['CUTS TIME'].append(float('NaN')) 
		elif len(data['CUTS TIME']) == it + 1:
			data['CUTS TIME'][-2] += data['CUTS TIME'][-1]
			data['CUTS TIME'] = data['CUTS TIME'][:-1]
		if len(data['CUTS CALLS']) < it:
			data['CUTS CALLS'].append(-1)
		elif len(data['CUTS CALLS']) == it + 1:
			data['CUTS CALLS'][-2] += data['CUTS CALLS'][-1]
			data['CUTS CALLS'] = data['CUTS CALLS'][:-1]
		if len(data['CUTS FOUND']) < it:
			data['CUTS FOUND'].append(-1)
		elif len(data['CUTS FOUND']) == it + 1:
			data['CUTS FOUND'][-2] += data['CUTS FOUND'][-1]
			data['CUTS FOUND'] = data['CUTS FOUND'][:-1]
		if len(data['CUTS APPLIED']) < it:
			data['CUTS APPLIED'].append(-1)
		elif len(data['CUTS APPLIED']) == it + 1:
			data['CUTS APPLIED'][-2] += data['CUTS APPLIED'][-1]
			data['CUTS APPLIED'] = data['CUTS APPLIED'][:-1]

		if len(data['FARKAS TIME']) < it:
			data['FARKAS TIME'].append(float('NaN'))
		if len(data['MASTER TIME']) < it:
			data['MASTER TIME'].append(float('NaN'))
		if len(data['PRICING TIME']) < it:
			data['PRICING TIME'].append(float('NaN'))

		if len(data['MIP OR KNAPSACK']) < it:
			data['MIP OR KNAPSACK'].append(-1)

		if len(data['DEGENERACY']) < it:
			data['DEGENERACY'].append(float('NaN')) 

		if len(data['CONS LINEAR']) < it:
			data['CONS LINEAR'].append(-1)
		if len(data['CONS KNAPSACK']) < it:
			data['CONS KNAPSACK'].append(-1)
		if len(data['CONS LOGICOR']) < it:
			data['CONS LOGICOR'].append(-1)
		if len(data['CONS SETPPC']) < it:
			data['CONS SETPPC'].append(-1)
		if len(data['CONS VARBOUND']) < it:
			data['CONS VARBOUND'].append(-1)

		if len(data['NBLOCKS']) < it:
			data['NBLOCKS'].append(-1)

		if len(data['AGGREGATIONS']) < it:
			data['AGGREGATIONS'].append(-1) 

		if len(data['SOLUTIONS FOUND']) < it:
			data['SOLUTIONS FOUND'].append(-1)
		if len(data['SOLUTIONS FIRST']) < it:
			data['SOLUTIONS FIRST'].append(float('NaN'))
		if len(data['SOLUTIONS BEST']) < it:
			data['SOLUTIONS BEST'].append(float('NaN'))
		if len(data['PD INTEGRAL']) < it:
			data['PD INTEGRAL'].append(float('NaN'))

		if len(data['MASTER NCONSS']) < it:
			data['MASTER NCONSS'].append(-1)
		if len(data['MASTER NVARS']) < it:
			data['MASTER NVARS'].append(-1)
		if len(data['LINKING VARS']) < it:
			data['LINKING VARS'].append(float('NaN')) 
		
		if len(data['BNB TREE NODES']) < it:
			data['BNB TREE NODES'].append(-1)
		if len(data['BNB TREE LEFT']) < it:
			data['BNB TREE LEFT'].append(-1)
		if len(data['BNB TREE DEPTH']) < it:
			data['BNB TREE DEPTH'].append(-1)

		if len(data['LP CALLS']) < it:
			data['LP CALLS'].append(-1)
		if len(data['LP TIME']) < it:
			data['LP TIME'].append(float('NaN'))
		if len(data['LP ITERATIONS']) < it:
			data['LP ITERATIONS'].append(-1)
		
		if len(data['LP FILE']) < it:
			data['LP FILE'].append(-1)
		if len(data['DEC FILE']) < it:
			data['DEC FILE'].append(-1)

#for key in data:
#	print("{0}: {1}".format(key, len(data[key])))

# build pandas data frame
pd.set_option("max_columns", 999)
df = pd.DataFrame(index=index, data=data)
df.to_pickle('pickles/' + outfile.split('/')[-1].replace('.out', '.pkl'))
#df.to_csv('pickles/' + outfile.replace('.out', '.csv'))

plot = False
if plot:
	s = "MASTER NVARS"
	plt.title(s)
	plt.xlabel('Instance')
	plt.ylabel(s)
	plt.xticks(range(len(index)), index, rotation='vertical')
	plt.tight_layout()
	plt.plot(df[s], 'x')
	plt.savefig('plots/{0}'.format(s.lower().replace(' ', '')))

twin = False 
if twin:
	s1 = "LP TIME"
	s2 = "PRICING TIME"
	
	# configuration
	fig, ax1 = plt.subplots(figsize=(15, 10))
	ax2 = ax1.twinx()
	ax1.set_xlabel("Instance")
	ax1.set_ylabel(s1, color='b') 
	ax2.set_ylabel(s2, color='r') 
	ax1.set_xticks(range(len(df[s1])))
	ax1.set_xticklabels(index, rotation=90)
	ax2.set_xticklabels(index, rotation=90)

	ax1.plot(df[s1].values, 'b.')
	ax2.plot(df[s2].values, 'r.')
	
	fig.tight_layout()
	plt.savefig('plots/{0}'.format(s1.lower().replace(' ', '') + s2.lower().replace(' ', '')))

# debugging output
def mok(k):
	if k == 0:
		return "none"
	elif k == 1:
		return "Knapsack"
	elif k == 2:
		return "MIP"
	elif k == 3:
		return "both"
	else:
		return "not given"

debug = False
if debug:
	print(df)
