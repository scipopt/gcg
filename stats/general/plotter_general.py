#!/usr/bin/env python3
import os
import sys
import re
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from sklearn.preprocessing import normalize


pickles = []
filters = []
I = []
totaltimes = []
numberofinstances = [1]

def read_pickle(pickle):
    pickles.append(pd.read_pickle(pickle))

def append_open_pickle(pickle):
    pickles.append(pickle)

def add_filter(attr, lb, ub):
    filters.append((attr, lb, ub))

def select_testset(name='filtered.test', decs=True, decs_strict=False):
    hd = open(name, 'w+')
    I = []
    for ts in range(len(pickles)):
        I.append([])
        pkl = pickles[ts]
        for index, row in pkl.iterrows():
            sat = True
            for filt in filters:
                if "%" in str(filt[1]) and "%" in str(filt[2]):
                    # percentage of TOTAL TIME. Only supported for times!
                    try:
                        lb = float(filt[1].replace("%",""))
                        ub = float(filt[2].replace("%",""))
                        percentage = float(row[filt[0]]) / float(row["TOTAL TIME"]) * 100
                        sat = sat and (lb <= percentage <= ub)
                    except:
                        # then TOTAL TIME was none, or a miscalculation occurred
                        sat = False
                else:
                    sat = sat and (filt[1] <= row[filt[0]] <= filt[2])
            if sat:
                if decs:
                    if row['LP FILE'] != -1 and row['DEC FILE'] != -1:
                        hd.write(row['LP FILE'] + ';' + row['DEC FILE'] + '\n')
                    elif row['LP FILE'] != -1:
                        hd.write(row['LP FILE'] + ';' + '\n')
                        print(f"Instance {index} (LP File {row['LP FILE']}) has no dec file")
                    else:
                        print("Instance {} has no lp or no dec file".format(index))
                else:
                    if row['LP FILE'] != -1:
                        hd.write(row['LP FILE'] + '\n')
                    else:
                        print("{} has no lp file".format(index))

def sum_pickles(names):
    print("Taking means...")
    npkls = len(pickles)
    for i in range(npkls):
        numberofinstances.append(len(pickles[i]["TOTAL TIME"]))
        totaltimes.append(pickles[i]["TOTAL TIME"].sum())
        pickles[i] = pickles[i].mean(axis=0,skipna=True,numeric_only=True).to_frame(name=str(names[i].split("/")[-1])).transpose()
    return npkls

def sum_mean_pickles(names):
    print("Taking means...")
    npkls = len(pickles)
    for i in range(npkls):
        numberofinstances.append(len(pickles[i]["TOTAL TIME"]))
        totaltimes.append(pickles[i]["TOTAL TIME"].sum())
        pickles[i] = pickles[i].sum(axis=0,skipna=True,numeric_only=True).to_frame(name=str(names[i].split("/")[-1])).transpose()
    return npkls

def pie(attributes, title='Plot', outdir="plots", filename="unknowntestset",m=[]):
    I = []
    ind = []
    A = {attr: [] for attr in attributes}
    for ts in range(len(pickles)):
        I.append([])
        pkl = pickles[ts]
        inst = 0
        for index, row in pkl.iterrows():
            sat = True
            for filt in filters:
                try:
                    sat = sat and (filt[1] <= row[filt[0]] <= filt[2])
                except:
                    print("Key '{}' not found. Please check spelling.".format(filt[0]))
                    exit()
            if sat:
                I[-1].append(inst)
                ind.append(index)
                for attr in attributes:
                    A[attr].append(row[attr])
            inst += 1

    labels = attributes.copy()
    for i in range(len(labels)):
        labels[i] = list((labels[i],0,np.inf))
        for filter in filters:
            #print("Comparing {} with {}".format(labels[i][0],filter[0]))
            if str(labels[i][0]) == str(filter[0]):
                #print("Found equality")
                if labels[i][1] < filter[1]:
                    #print("Setting {} to value {}".format(labels[i][1],filter[1]))
                    labels[i][1] = filter[1]
                if labels[i][2] > filter[2]:
                    #print("Setting")
                    labels[i][2] = filter[2]
        labels[i] = str(labels[i][0]) + " [" + str(labels[i][1]) + ", " + str(labels[i][2]) + "]"
        i += 1

    plt.figure(figsize=(15, 10))
    plt.title(str(title) + " of the testset " + "'" + str(filename) + "' with " + str(len(ind)) + " instances")
    #plt.xlabel('Instance')
    #plt.ylabel(', '.join(attributes))
    plt.xticks(range(len(ind)), ind, rotation='vertical')
    plt.tight_layout()
    i = 0
    msize = '8'
    if len(m) < len(attributes):
        if len(m) == 0 and len(attributes) == 2:
            m = ["_","|"]
            msize = '10'
        elif len(m) == 0 and len(attributes) == 3:
            m = [4,5,6]
            msize = '8'
        elif len(m) == 0 and len(attributes) == 4:
            m = [4,5,6,7]
            msize = '8'
        else:
            for j in range(len(attributes)-len(m)):
                m.append('x')
    print([A[a][0] for a in attributes])
    plt.pie([A[a][0] for a in attributes], labels=labels)

    plt.legend(loc='best',prop={'size': 12})

    plt.savefig(os.path.join(outdir,'{}.pie.pdf'.format(filename)))
    print("Saved graph '\033[1m\033[92m{}\033[0m' as '{}.pie.pdf'".format(title,filename))

def plot(attributes, title='Plot', outdir="plots", filename="unknowntestset",m=[]):
    I = []
    ind = []
    A = {attr: [] for attr in attributes}
    for ts in range(len(pickles)):
        I.append([])
        pkl = pickles[ts]
        inst = 0
        for index, row in pkl.iterrows():
            sat = True
            for filt in filters:
                try:
                    sat = sat and (filt[1] <= row[filt[0]] <= filt[2])
                except:
                    print("Key '{}' not found. Please check spelling.".format(filt[0]))
                    exit()
            if sat:
                I[-1].append(inst)
                ind.append(index)
                for attr in attributes:
                    A[attr].append(row[attr])
            inst += 1

    labels = attributes.copy()
    for i in range(len(labels)):
        labels[i] = list((labels[i],0,np.inf))
        for filter in filters:
            #print("Comparing {} with {}".format(labels[i][0],filter[0]))
            if str(labels[i][0]) == str(filter[0]):
                #print("Found equality")
                if labels[i][1] < filter[1]:
                    #print("Setting {} to value {}".format(labels[i][1],filter[1]))
                    labels[i][1] = filter[1]
                if labels[i][2] > filter[2]:
                    #print("Setting")
                    labels[i][2] = filter[2]
        labels[i] = str(labels[i][0]) + " [" + str(labels[i][1]) + ", " + str(labels[i][2]) + "]"
        i += 1

    plt.figure(figsize=(15, 10))
    plt.title(str(title) + " of the testset " + "'" + str(filename) + "' with " + str(len(ind)) + " instances")
    plt.xlabel('Instance')
    plt.ylabel(', '.join(attributes))
    plt.xticks(range(len(ind)), ind, rotation='vertical')
    plt.tight_layout()
    i = 0
    msize = '8'
    if len(m) < len(attributes):
        if len(m) == 0 and len(attributes) == 2:
            m = ["_","|"]
            msize = '10'
        elif len(m) == 0 and len(attributes) == 3:
            m = [4,5,6]
            msize = '8'
        elif len(m) == 0 and len(attributes) == 4:
            m = [4,5,6,7]
            msize = '8'
        else:
            for j in range(len(attributes)-len(m)):
                m.append('x')
    for attr in attributes:
        plt.plot(A[attr], label=labels[i], marker=m[i],linestyle='',markersize=msize)
        i += 1
    plt.legend(loc='best',prop={'size': 12})

    plt.savefig(os.path.join(outdir,'{}.plot.pdf'.format(filename)))
    print("Saved graph '\033[1m\033[92m{}\033[0m' as '{}.plot.pdf'".format(title,filename))

def twinplot(attributes1, attributes2, title='Twinplot', outdir="plots", filename="unknowntestset"):
    I = []
    ind = []
    A1 = {attr: [] for attr in attributes1}
    A2 = {attr: [] for attr in attributes2}
    for ts in range(len(pickles)):
        I.append([])
        pkl = pickles[ts]
        inst = 0
        for index, row in pkl.iterrows():
            sat = True
            for filt in filters:
                try:
                    sat = sat and (filt[1] <= row[filt[0]] <= filt[2])
                except:
                    print("Key '{}' not found. Please check spelling.".format(filt[0]))
                    exit()
            if sat:
                I[-1].append(inst)
                ind.append(index)
                for attr in attributes1:
                    A1[attr].append(row[attr])
                for attr in attributes2:
                    A2[attr].append(row[attr])
            inst += 1

    labels1 = attributes2.copy()
    for i in range(len(labels1)):
        labels1[i] = list((labels1[i],0,np.inf))
        for filter in filters:
            #print("Comparing {} with {}".format(labels1[i][0],filter[0]))
            if str(labels1[i][0]) == str(filter[0]):
                #print("Found equality")
                if labels1[i][1] < filter[1]:
                    #print("Setting {} to value {}".format(labels1[i][1],filter[1]))
                    labels1[i][1] = filter[1]
                if labels1[i][2] > filter[2]:
                    #print("Setting")
                    labels1[i][2] = filter[2]
        labels1[i] = str(labels1[i][0]) + " [" + str(labels1[i][1]) + ", " + str(labels1[i][2]) + "]"
        i += 1

    labels2 = attributes1.copy()
    for i in range(len(labels2)):
        labels2[i] = list((labels2[i],0,np.inf))
        for filter in filters:
            #print("Comparing {} with {}".format(labels2[i][0],filter[0]))
            if str(labels2[i][0]) == str(filter[0]):
                #print("Found equality")
                if labels2[i][1] < filter[1]:
                    #print("Setting {} to value {}".format(labels2[i][1],filter[1]))
                    labels2[i][1] = filter[1]
                if labels2[i][2] > filter[2]:
                    #print("Setting")
                    labels2[i][2] = filter[2]
        labels2[i] = str(labels2[i][0]) + " [" + str(labels2[i][1]) + ", " + str(labels2[i][2]) + "]"
        i += 1

    fig, ax1 = plt.subplots(figsize=(15, 10))
    plt.title(str(title) + " of the testset " + "'" + str(filename) + "' with " + str(len(ind)) + " instances")
    ax2 = ax1.twinx()
    ax1.set_xlabel("Instance")
    ax1.set_ylabel(", ".join(attributes1))
    ax2.set_ylabel(", ".join(attributes2))
    ax1.set_xticks(range(len(ind)))
    ax1.set_xticklabels(ind, rotation='vertical')
    fig.tight_layout()
    i = 0
    for attr in attributes1:
        ax1.plot(A1[attr], label=labels1[0], marker='x',linestyle='',color='r')
        i += 1
    for attr in attributes2:
        ax2.plot(A2[attr], label=labels2[0], marker='x',linestyle='',color='b')
        i += 1

    ax1.legend(loc='upper left',prop={'size': 12})
    ax2.legend(loc='upper right',prop={'size': 12})

    plt.savefig(os.path.join(outdir,'{}.twin.pdf'.format(filename)))
    print("Saved graph '\033[1m\033[92m{}\033[0m' as '{}.twin.pdf'".format(title,filename))

def bubbleplot(attribute1, attribute2, title='Bubbleplot', outdir="plots", filename="unknowntestset", proximity=[10,10]):
    I = []
    ind = []
    A = []
    B = []
    for ts in range(len(pickles)):
        I.append([])
        pkl = pickles[ts]
        inst = 0
        for index, row in pkl.iterrows():
            sat = True
            for filt in filters:
                try:
                    sat = sat and (filt[1] <= row[filt[0]] <= filt[2])
                except:
                    print("Key '{}' not found. Please check spelling.".format(filt[0]))
                    exit()
            if sat:
                I[-1].append(inst)
                ind.append(index)
                A.append(row[attribute1])
                B.append(row[attribute2])
            inst += 1

    labels = []
    labels.append(list((attribute1,0,np.inf)))
    labels.append(list((attribute2,0,np.inf)))
    for i in range(len(labels)):
        for filter in filters:
            #print("Comparing {} with {}".format(labels[i][0],filter[0]))
            if str(labels[i][0]) == str(filter[0]):
                #print("Found equality")
                if labels[i][1] < filter[1]:
                    #print("Setting {} to value {}".format(labels[i][1],filter[1]))
                    labels[i][1] = filter[1]
                if labels[i][2] > filter[2]:
                    #print("Setting")
                    labels[i][2] = filter[2]
        labels[i] = str(labels[i][0]) + " [" + str(labels[i][1]) + ", " + str(labels[i][2]) + "]"
        i += 1


    # if points are close up to x_prox/y_prox units, they are cumulated
    prox = [1, 1]
    if len(proximity) == 1:
        prox[0] = proximity[0]
        prox[1] = proximity[0]
    elif len(proximity) == 2:
        prox = proximity
    else:
        print("Given proximity was erroneous.")
        prox[0] = 10
        prox[1] = 10
    print("Using proximity {}".format(prox))

    # determine bubble positions and sizes
    AA = []
    BB = []
    CC = []
    for i in range(len(A)):
        ain = False
        for j in range(len(AA)):
            if abs(A[i]-AA[j]) <= prox[0] and abs(B[i]-BB[j]) <= prox[1]:
                #print("Found close points!")
                CC[j] += 1
                ain = True
        if not ain:
            AA.append(A[i])
            BB.append(B[i])
            CC.append(1)

    # increase bubble size a bit for all bubbles
    for i in range(len(CC)):
        CC[i] = 4*CC[i]**2

    plt.figure(figsize=(10, 10))
    plt.figtext(.01,.01,'Proximity factors: '+str(prox))
    plt.title(str(title) + " of the testset " + "'" + str(filename) + "' with " + str(len(ind)) + " instances")
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.tight_layout()
    plt.scatter(AA, BB, s=CC, alpha=.5)

    plt.savefig(os.path.join(outdir,'{}.bubble.pdf'.format(filename)))
    print("Saved graph '\033[1m\033[92m{}\033[0m' as '{}.bubble.pdf'".format(title,filename))

def time_distribution(times, nbuckets, title='Time Distribution', type='bar', outdir="plots", filename="unknowntestset", linearbucketing=False):
    I = []
    ind = []
    singleinstancename=""
    #print("   Using arguments {}".format(times))
    # use same colors as in GCG decomp visualizations
    COLOR_WHITE = "#FFFFFF"
    COLOR_BLUE1 = "#ACBCE9"
    COLOR_BLUE2 = "#718CDB"
    COLOR_BLUE3 = "#3C64DD"
    COLOR_BLUE4 = "#1340C7"
    COLOR_BLUE5 = "#1F377D"
    COLOR_ORANGE1 = "#FFD88F"
    COLOR_ORANGE2 = "#FFCB69"
    COLOR_ORANGE3 = "#FFB72D"
    COLOR_BROWN1 = "#B38208"
    COLOR_BROWN2 = "#886100"
    COLOR_BROWN3 = "#443000"
    COLOR_BLACK = "#000000"

    colormapping = [
    ("TOTAL TIME",    COLOR_BLACK),
    ("DETECTION TIME",  COLOR_BLUE3),
    ("RMP LP TIME",     COLOR_BLUE2),
    ("PRICING TIME",    COLOR_ORANGE3),
    ("PRICING SOLVER TIME",    COLOR_ORANGE2),
    ("CUTS TIME ORIG",       COLOR_ORANGE2),
    ("CUTS TIME MASTER",       COLOR_ORANGE2),
    ("HEURISTICS TIME ORIG", COLOR_ORANGE1),
    ("HEURISTICS TIME MASTER", COLOR_ORANGE1),
    ("FARKAS TIME",     COLOR_BROWN3),
    ("ROOT NODE TIME",  COLOR_BROWN2),
    ("ORIGINAL LP TIME",COLOR_BROWN1),
    ("COPYING TIME",    '0.1'),
    ("READING TIME",    '0.2'),
    ("PRESOLVING TIME", '0.3'),

    ("BRANCHING RULE CALLS GENERIC",  COLOR_BLUE2),
    ("BRANCHING RULE CALLS ORIG",  COLOR_BROWN2),
    ("BRANCHING RULE CALLS RELPSPROB",     COLOR_ORANGE2),
    ("BRANCHING RULE CALLS RYANFOSTER",       COLOR_BLACK),
    ]

    colors = []

    # Set color mapping
    for time in times:
        # If the time was defined in the colormapping var, assign this color
        if time in (np.transpose(colormapping))[0]:
            colors.append(colormapping[list(np.transpose(colormapping)[0]).index(time)][1])
        # Else, user probably wants to highlight her/his specific time, use red
        else:
            colors.append('#ff0000')

    A = []
    for i in range(len(times)):
        A.append([])

    for ts in range(len(pickles)):
        I.append([])
        pkl = pickles[ts]
        inst = 0
        for index, row in pkl.iterrows():
            sat = True
            for filt in filters:
                try:
                    sat = sat and (filt[1] <= row[filt[0]] <= filt[2])
                    #if not sat: print("{} does not satisfy filter {} with its value {}.".format(index,filt,row[filt[0]]))
                except:
                    print("Key '{}' not found. Please check spelling.".format(filt[0]))
                    exit()
            if sat:
                I[-1].append(inst)
                ind.append(index)
                for i in range(len(times)): # insert the times for one instance into the array
                    #print("Round {}".format(i))
                    A[i].append(row[times[i]])
            inst += 1

    A = np.array(A)

    if type == "compare":
        settings_place=-5
    else:
        settings_place=-3

    if ind == [] and not type == "compare":
        print("Empty pickle file or no instances left after filtering.\nTerminating.")
        exit()
    else:
        ind_temp = ind.copy()
        ind_with_newlines = ind.copy()
        if ind[0].startswith("check."):
            # Then its a testset and we need the second place
            singleinstancename = ind[0].split(".")[1]+ "." + ind[0].split(".")[settings_place]
            for i in range(len(ind)):
                ind_temp[i] = ind[i].split(".")[1] + "." + ind[i].split(".")[settings_place]
                ind_with_newlines[i] = ind[i].split(".")[1] + "\n" + ind[i].split(".")[settings_place] + "\n(" + str(numberofinstances[i+1]) + " instances)"
            ind = ind_temp
        else:
            print("You modified your out file names. Unexpected behavior might occur.")



    # FARKAS is included in PRICING
    if "FARKAS TIME" in times and "PRICING TIME" in times:
        for n in range(len(A[0])):
            A[times.index("PRICING TIME")][n] -= A[times.index("FARKAS TIME")][n]

    if "PRICING SOLVER TIME" in times and "PRICING TIME" in times:
        for n in range(len(A[0])):
            A[times.index("PRICING TIME")][n] -= A[times.index("PRICING SOLVER TIME")][n]

    # READING/COPYING/DETECTING is not included in OVERALL
    if "READING TIME" in times:
        for n in range(len(A[0])):
            A[times.index("TOTAL TIME")][n] += A[times.index("READING TIME")][n]
    if "COPYING TIME" in times:
        for n in range(len(A[0])):
            A[times.index("TOTAL TIME")][n] += A[times.index("COPYING TIME")][n]
    if "DETECTION TIME" in times:
        for n in range(len(A[0])):
            A[times.index("TOTAL TIME")][n] += A[times.index("DETECTION TIME")][n]
    if "ORIGINAL LP TIME" in times:
        for n in range(len(A[0])):
            A[times.index("TOTAL TIME")][n] += A[times.index("ORIGINAL LP TIME")][n]


    # Warn if MASTER is included
    if "MASTER TIME" in times:
        # If MASTER TIME should be plotted at some point, subtract PRICING, both LP, HEURISTICS and CUTS times from it
        print("    Error: MASTER TIME should not be plotted. Use RMP LP TIME instead.\n    Terminating.")
        exit()

    if "ROOT NODE TIME" in times:
        if "PRICING TIME" in times:
            print("    Error: ROOT NODE TIME cannot be plotted with PRICING TIME, since it is partly contained in PRICING.\n    Terminating.")
            exit()
        if "RMP LP TIME" in times:
            print("    Error: ROOT NODE TIME cannot be plotted with RMP LP TIME, since it is partly contained in PRICING.\n    Terminating.")
            exit()

    # Check, if overall time >= sum of all other times
    n = 0
    terrorn = 0
    terror = False
    while n < len(A[0])-1:
        if (A[0][n] - sum(A[i][n] for i in range(1,len(A))) < 0):
            print("    Warning: Time measurement error found for problem {}, set {}. Deleted entry.".format(n, [row[n] for row in A]))
            A[0][n] = 0
            A = np.delete(A,n,1)
            terror = True
            terrorn += 1
        else:
            n += 1
    if terror:
        if "ROOT NODE TIME" in times:
            print("    Information: ROOT NODE TIME should not be plotted with any\
             other time except for READING, COPYING and DETECTION.")
        print("    {} time measurement errors found. Plot will have no expressive power.\n    Terminating.".format(terrorn))
        exit()
    xlabel = None

    old_ind = ind.copy()
    # Grouping of bars before processing
    if type == "grouped_bar":
        # Bucket linearly (don't care how many instances per group, just linspace it)
        if linearbucketing:
            print("    Information: Linear bucketing is enabled.")
            if nbuckets == -1:
                nbuckets = np.sqrt((np.ceil(np.amax(A,axis=1)[0]) - np.floor(np.amin(A,axis=1)[0])))
                if nbuckets < 20:
                    if len(A[0]) > 19:
                        nbuckets = 20
                    else:
                        nbuckets = len(A[0])
            # Define Buckets
            buckets = np.linspace(np.floor(np.amin(A,axis=1)[0]),np.ceil(np.amax(A,axis=1)[0]),nbuckets+1)
            xlabel = "CPU time in steps of " + "{0:.0f}".format(round(buckets[1] - buckets[0]),2) + " (" + str(len(buckets)/10) + "-percentiles)"
        # Bucket dependent of #instances. Each bucket will get a more or less equal amount of instances
        # Soft-limit, since computation only with 2 decimals, SCIP calculates with 6
        # Will need more time than linear bucketing
        else:
            # first sort array
            A = zip(*A)
            A = sorted(A,key=itemgetter(0)) # 0 for sort by other, 1 for sort by master, ...
            A = list(zip(*A))
            # initialize buckets
            buckets = []
            # If not set, define a sensible amount of buckets
            if nbuckets == -1:
                if len(A[0]) < 10:
                    nbuckets = len(A[0])
                else:
                    nbuckets = 10
            else:
                nbuckets+=0

            instances_per_bucket = int(np.floor(len(A[0]) / nbuckets))
            if instances_per_bucket == 0:
                print("More buckets than instances chosen. Using default value.")
                if len(A[0]) < 10:
                    nbuckets = len(A[0])
                    instances_per_bucket = int(np.floor(len(A[0]) / nbuckets))
                else:
                    nbuckets = 10
                    instances_per_bucket = int(np.floor(len(A[0]) / nbuckets))
            # loop and take every instance_per_bucket'th element
            n = 0
            while n < len(A[0]) and len(buckets) < nbuckets:
                #print(n)
                #print("mod")
                #print(instances_per_bucket)
                if (n % instances_per_bucket) == 0:
                    buckets.append(np.around(A[0][n],decimals=2))
                n += 1
            buckets.append(np.around(A[0][len(A[0])-1],decimals=2))

        # Check if Bucketing works
        #print("   Bars grouped into {}.".format(buckets))

        # get the time of the last tick before normalization
        if not linearbucketing:
            endtick = [np.around(A[0][len(A[0])-1],decimals=2)]

        # Put instances into buckets
        Bdic = {lowerBound: [] for lowerBound in buckets}
        A=list(zip(*A))
        for i in range(len(A)):
            for j in range(len(buckets))[:len(buckets)-1]:
                try:
                    if (A[i][0] <= buckets[j+1]) and (A[i][0] >= buckets[j]):
                        Bdic[buckets[j]].append(A[i])
                        #print("Inserted {} into Bucket {}".format(A[i][:], buckets[j]))
                        break
                except TypeError:
                    print("TypeError")
        B = np.array([[]])

        # See if the grouping was roughly equal (last one always zero, since it is only the upper bound, it represents no bucket)
        print("   with bucket sizes {}".format([len(Bdic[key]) for key in Bdic][:-1]))

        # Take means of each bucket, if bucket empty delete it
        for key in list(Bdic.keys()):
            if not Bdic[key] == []:
                B=np.append(B,np.mean(np.array(Bdic[key]), axis=0))
            else:
                del Bdic[key]
        try:
            A = np.reshape(B,(len(Bdic),len(times)))
        except:
            print("BucketException")
            exit()
        A = np.array(list(zip(*A)))
        ind = np.round(sorted(list(Bdic.keys())),decimals=2)
        if not linearbucketing:
            xlabel = "CPU time (~{0:.0f} instances per bar)".format(round(sum([len(Bdic[key]) for key in Bdic][:-1])/len([len(Bdic[key]) for key in Bdic][:-1]))+1)

    # For the comparison pie plot, we want to add a label with the total time,
    # so we copy it here before normalization
    averagetimes = [sum(A[i]) for i in range(len(A))]

    # Subtract every other item from the overall time to get "OTHER"
    for n in range(len(A[0])):
        A[0][n] = A[0][n] - sum(A[i][n] for i in range(1,len(times)))
        if A[0][n] < 0:
            print("Measurement/time inclusion error possible. Treat results with caution.")
            A[0][n] = 0

    # Make matrix stochastic and sort it by the first argument (default: PRICING TIME)
    A = normalize(A, axis=0, norm='l1') # Axis must be 1

    if not type == "grouped_bar":
        # Cast as object, to add string column
        A = np.array(A,dtype=object)
        # Append instance names (string column)
        A = np.append(A,[ind], axis=0)

    if type == "bar" or type == "grouped_bar":# Sort array (with instance names)
        A = np.array(list(zip(*A)),dtype=object)
        A = np.array(sorted(A,key=itemgetter(1))) # 0 for sort by other, 1 for sort by first argument in -t, ...
        A = np.array(list(zip(*A)))

    if len(ind) == 1 or (len(ind) < 10 and (type == "compare" or type == "comparepie")): singleinstance = True
    else: singleinstance = False

    if not type == "grouped_bar":
        # Get instance names back
        ind = A[:][len(times)]
        # Remove instance names
        A = np.delete(A,len(times),axis=0).astype(np.float)

    # Prepare the plot
    if singleinstance:
        plt.figure(figsize=(15, 10))
        if type == "pie":
            plt.title(str(filename) + ' (' + str(numberofinstances[0]) + " instances)")
        else:
            plt.title(str(singleinstancename) + ' (' + str(numberofinstances[0]) + " instances)")
        if type == "comparepie" or type == "compare":
            plt.title("Comparison")
    else:
        plt.figure(figsize=(25, 10))
        plt.title(title + " of the testset " + "'" + str(filename) + "' with " + str(len(old_ind)) + " instances")
    if not (type == "pie" or type == "comparepie"):
        plt.ylabel(times[0])

    # Label corresponding to bucketed/not bucketed
    if not xlabel and not (type == "compare" or type == "comparepie" or type == "pie" or (type == "bar" and len(ind) == 1)):
        plt.xlabel("Instances")
    else:
        plt.xlabel(xlabel)

    # only plot instance names if they can be read
    if (type == "grouped_bar") and not linearbucketing:
        ind = list(ind)+endtick
        plt.xticks(range(len(ind)), ind, rotation='horizontal')
    elif type == "compare":
        plt.xticks(range(len(ind)), ind_with_newlines)
    elif len(ind) < 125:
        if len(ind) < 5:
            plt.xticks(range(len(ind)), ind, rotation='horizontal')
        else:
            plt.xticks(range(len(ind)), ind, rotation='vertical')

    ind_markers = ind.copy()
    ind = np.arange(len(A[0]))
    A = np.array(A)

    labels = times.copy()
    for i in range(len(labels)):
        labels[i] = list((labels[i],0,np.inf))
        for filter in filters:
            #print("Comparing {} with {}".format(labels[i][0],filter[0]))
            if str(labels[i][0]) == str(filter[0]):
                #print("Found equality")
                if labels[i][1] < filter[1]:
                    #print("Setting {} to value {}".format(labels[i][1],filter[1]))
                    labels[i][1] = filter[1]
                if labels[i][2] > filter[2]:
                    #print("Setting")
                    labels[i][2] = filter[2]
        labels[i] = str([str(labels[i][0]) if labels[i][0] != "TOTAL TIME" else "OTHER"][0]) + " [" + str(labels[i][1]) + ", " + str(labels[i][2]) + "]"
        i += 1

    # Plotting a bargraph
    if type == "grouped_bar":
        # Plot all times. For color scheme: check if color is RGBA (i.e. a tuple) or just a string/number
        for i in range(1,len(times)):
            plt.bar(ind+0.5, A[i], bottom=A[:][1:i].sum(axis=0) if i > 1 else 0, label=labels[i], color=[eval(colors[i]) if colors[i][0] == "(" else colors[i]])
        plt.bar(ind+0.5, A[0], bottom=A[:][1:len(times)].sum(axis=0), label=labels[0], color=colors[0]) # Plot the rest until 1

    if type == "bar":
        # Plot all times
        for i in range(1,len(times)):
            plt.bar(ind, A[i], bottom=A[:][1:i].sum(axis=0) if i > 1 else 0, label=labels[i], color=[eval(colors[i]) if colors[i][0] == "(" else colors[i]])
        plt.bar(ind, A[0], bottom=A[:][1:len(times)].sum(axis=0), label=labels[0], color=colors[0]) # Plot the rest until 1

    if type == "compare":
        # Plot all times
        for i in range(1,len(times)):
            plt.bar(ind, A[i], bottom=A[:][1:i].sum(axis=0) if i > 1 else 0, label=labels[i], color=[eval(colors[i]) if colors[i][0] == "(" else colors[i]])
        plt.bar(ind, A[0], bottom=A[:][1:len(times)].sum(axis=0), label=labels[0], color=colors[0]) # Plot the rest until 1

    # Plotting a simple plot
    elif type == "plot":
        plt.plot(ind, np.ones(len(A[0])))
        for i in range(1,len(times)):
            plt.plot(ind, A[i], 'x', label=labels[i], color=eval(colors[i]) if colors[i][0] == "(" else colors[i])
        plt.plot(ind, A[0], 'x', label=labels[0], color=colors[0])

    elif type == "pie":
        if len(A[0]) > 1:
            print("More than one instance cannot be plotted with a simple pie chart.\nTerminating.")
            exit()
        _, _, autotexts = plt.pie(A.flatten(), labels=labels, colors=[(eval(colors[i]) if colors[i][0] == "(" else colors[i]) for i in range(len(colors))], autopct='%1.1f%%', radius=0.8, textprops=dict(size="large"),shadow=True)
        for autotext in autotexts:
            autotext.set_color('white')#plt.pie(ind, A[0], labels="OTHER", colors=colors[0]) # Plot the rest until 1
        plt.figtext(.01,.01,'The total runtime of the testset (' + str(numberofinstances[0]) + ' instances) was ' + str("{0:.2f}".format(totaltimes[0])) + 's.',size='medium')
        plt.figtext(.01,.01+1/60,'The average runtime of an instance was ' + str("{0:.2f}".format(averagetimes[0])) + 's.',size='medium')


    elif type == "comparepie":
        A=A.transpose()
        size=1/(len(A)+1)
        for i in range(len(A)):
            plt.pie(A[i], labels=[labels[i] if labels[i] != "TOTAL TIME" else "OTHER" for i in range(len(labels))] if i == 0 else None, colors=[(eval(colors[i]) if colors[i][0] == "(" else colors[i]) for i in range(len(colors))], radius=1-size*i,wedgeprops=dict(width=size, edgecolor='w'),startangle=90,shadow=True)
            plt.figtext(.01,.01+i/60,'The total runtime of testset ' + ind_markers[i] + ' (#' + str(i+1) + ' from outside, ' + str(numberofinstances[i+1]) + ' instances) was ' + str("{0:.2f}".format(totaltimes[i])) + 's.',size='medium')
        #plt.pie(ind, A[0], labels="OTHER", colors=colors[0]) # Plot the rest until 1

    # Set legend, layout
    if type == "comparepie" or type == "pie":
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.gca().legend(reversed(handles), reversed(labels), loc='upper left')
    else:
        plt.legend(loc='best')

    plt.tight_layout()
    # remove settings for filename if compare(pie)
    if type == "compare" or type == "comparepie":
        filename = str(filename).split('.')[0]
        singleinstancename = str(singleinstancename).split('.')[0]

    # Save the plots
    if not singleinstance:
        plt.savefig(os.path.join(outdir,'{}.timedist.{}.pdf'.format(filename,type)))
        print("Saved graph '\033[1m\033[92m{}\033[0m' as '{}.timedist.{}.pdf'".format(type,filename,type))
    else:
        if type == "bar": type = "singlebar"
        plt.savefig(os.path.join(outdir,'{}.timedist.{}.pdf'.format(singleinstancename,type)))
        print("Saved graph '\033[1m\033[92m{}\033[0m' as '{}.timedist.{}.pdf'".format(type,singleinstancename,type))
    #plt.show()

if __name__ == '__main__':
    print("Warning: This script should not be executed directly. Please use plot.py, bubble.py, time.py or twin.py.\nTerminating.")
    exit()
