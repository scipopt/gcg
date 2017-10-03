#!/usr/bin/env python2

import sys
import os
import argparse

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

# Define global variables
params = {}

def parse_arguments(args):
    """
    Parse the command-line arguments
    :param args: Command-line arguments provided during execution
    :return: Parsed argumetns
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-n', '--allnodes', action='store_true',
                        help='set this flag, to plot the pricing-statistics for all nodes, not just the root')

    parser.add_argument('-s', '--splitrounds', type=int,
                        default=0,
                        help='the plot can be split into pieces containing a maximum of SPLITROUNDS rounds each (default is no splitting)')

    parser.add_argument('filename', nargs='+',
                        help='Name of the files to be used for the bound plots')

    parsed_args = parser.parse_args(args)
    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['outdir'] = args.outdir
    params['root_only'] = not args.allnodes
    params['nSplit'] = args.splitrounds

def get_colmap(pricers):
    """
    Returns a list of colors, with same length as pricers, that can be used for the bar plot
    Each pricing problem has its own color
    :param pricers: a list with pricing_problem ids
    :return: a list of colors as demanded by pyplot.bar()
    """
    pricer_to_color = {}
    col_ind = 0

    # build the mapping pricer id -> color id
    for p in pricers:
        if not p in pricer_to_color:
            pricer_to_color[p] = col_ind
            col_ind += 1

    # get a color map of the right length, so that each color-id gets its own color
    cmap = plt.get_cmap('nipy_spectral',len(pricer_to_color))

    # build a list of colors and return it
    colors = [cmap(pricer_to_color[p]) for p in pricers]
    return colors

def make_plots(data, name):
    """
    Make the plots from the structured data
    :param data: dataframe with the collected data
    :param name: the problemname
    :return:
    """
    # set parameters for the plot
    ymin = -0.15

    # flat out the data again (maybe in the final version, there should be no dataframe to begin with)
    flat_data = data.reset_index()
    flat_data.time = flat_data.time + 0.01 # workaround for most times being zero (0,01s is SCIPs smallest time-interval)

    # define position, width and height of the peaks; the former are defined by time, the latter by nVars
    x = (flat_data.time.cumsum()-flat_data.time).values
    y = (flat_data.nVars - ymin).values
    widths = flat_data.time.values
    colors = get_colmap(flat_data['pricing_prob'].values)

    # make the bar plot
    plt.bar(x, y, widths, bottom = ymin, align = 'edge', edgecolor = 'k', color = colors, label='pricing problems')

    # formatting
    ax = plt.gca()
    ymax = max(y)
    ax.set_xlim([0,x[-1]+widths[-1]])
    ax.set_ylim([ymin,ymax])
    old_yticks = ax.get_yticks()
    new_yticks = []
    for i,n in enumerate(old_yticks):
        if abs(n - int(n)) < 0.01 and n >= 0:
            new_yticks.append(n)
    ax.set_yticks(new_yticks)
    ax.set_xlabel('time')
    ax.set_ylabel('nVars')
    plt.gcf().subplots_adjust(top=0.835)

    # add information about the pricing rounds
    prev_rnd = min(flat_data['pricing_round'])
    prev_x = 0
    texts = []
    for i in range(len(x)):
        rnd = flat_data['pricing_round'][i]
        if rnd > prev_rnd:
            ax.plot([x[i],x[i]],[ymin,ymax],'r--',linewidth=0.8)
            texts.append(ax.text((x[i] + prev_x)/2. + 0.0033, ymax*1.01, 'Round '+str(prev_rnd), rotation='vertical',va='bottom', ha='center'))
            prev_rnd = rnd
            prev_x = x[i]
    texts.append(ax.text((x[-1] + widths[-1] + prev_x)/2. + 0.0033, ymax*1.01, 'Round '+str(prev_rnd), rotation='vertical',va='bottom', ha='center'))

    # check for overlapping texts
    rend = plt.gcf().canvas.get_renderer()
    for i,txt in enumerate(texts):
        if not txt.get_visible():
            continue
        bb = txt.get_window_extent(renderer=rend)
        for nxt_txt in texts[(i+1):]:
            nxt_bb = nxt_txt.get_window_extent(renderer=rend)
            if bb.overlaps(nxt_bb):
                nxt_txt.set_visible(False)
            else:
                break

    # save the figure
    plt.savefig(params['outdir'] + '/' + name + '.png')
    plt.close()

    print '    saved figure'

def generate_files(files):
    """
    Parse the files and structure the pricing-data in a dataframe then make the plots
    :param files: List of files to be parsed
    :return:
    """
    for file in files:
        with open(file) as _file:
            # initialize all index-lists
            ind_node = []
            ind_pricing_round = []
            ind_stab_round = []
            ind_pricing_prob = []

            # initialize the value-lists
            val_time = []
            val_nVars = []

            # initialize all counters
            node = 0
            pricing_round = 0
            stab_round = 0
            pricing_prob = 0

            # initialize all other variables
            problemFileName = None
            done = False

            for line in _file:
                # if the file is a out-file, generated by the check-script, reset the variables whenever a new instance starts
                if line.startswith("@01"):
                    if problemFileName and not done:
                        print '    could not parse the data...'

                    # initialize all index-lists
                    ind_node = []
                    ind_pricing_round = []
                    ind_stab_round = []
                    ind_pricing_prob = []

                    # initialize the value-lists
                    val_time = []
                    val_nVars = []

                    # initialize all counters
                    node = 0
                    pricing_round = 0
                    stab_round = 0
                    pricing_prob = 0

                    # initialize all other variables
                    problemFileName = None
                    done = False

                # if the problem is already processed, continue
                elif done:
                    continue

                elif not problemFileName and line.startswith("read problem "):
                    # get the problem name from the file name as in "check.awk"
                    tmparray = line.split("<")[-1].replace(">","").replace("\n","").split("/")[-1].split(".")
                    problemFileName = tmparray[0]
                    if tmparray[-1] == "gz" or tmparray[-1] == "z" or tmparray[-1] == "GZ" or tmparray[-1] == "Z":
                        tmparray.pop()
                    for i in range(1,len(tmparray)-1):
                        problemFileName += "." + tmparray[i]
                    print 'entering', problemFileName

                # pricer statistics end
                elif line.startswith("SCIP Status        :"):
                    # continue if no data is found
                    if not ind_node or not ind_pricing_round or not ind_stab_round or not ind_pricing_prob or not val_time or not val_nVars:
                        print '    no pricing data found'
                        done = True
                        continue

                    index = pd.MultiIndex.from_arrays([ind_node, ind_pricing_round, ind_stab_round, ind_pricing_prob],
                                                      names=["node", "pricing_round", "stab_round", "pricing_prob"])
                    data = {'time': val_time, 'nVars': val_nVars}
                    df = pd.DataFrame(data=data, index = index)

                    # split the data into pieces of params['nSplit'] rounds
                    # do nothing if the paramter nSplit is not set (equals zero)
                    if params['nSplit'] <= 0:
                        make_plots(df, problemFileName)
                    else:
                        maxRnd  = max(df.index.get_level_values('pricing_round').unique().values)
                        fromRnd = 0
                        for i in range(1,maxRnd+1):
                            if i % params['nSplit'] <> 0 and i <> maxRnd:
                                continue
                            toRnd = i
                            make_plots(df.query('@fromRnd < pricing_round <= @toRnd'), problemFileName + '_rounds' + str(fromRnd + 1) + 'to' + str(toRnd))
                            fromRnd = toRnd

                    done = True

                    print '    leaving', problemFileName
                    continue

                # ignore all other lines, that do not contain pricer statistics messages
                elif not line.startswith("[src/pricer_gcg.cpp:"):
                    continue

                # ignore lines, where the output ends abrubtly (e.g. when the hard limit of the check-script is reached)
                if line.find("@0")<>-1:
                    continue

                # extract the pricing-statistics message
                message = line.split("] statistic: ")[-1]

                if message.startswith("New pricing round at node"):
                    node = int(message.split()[-1])
                    pricing_round += 1

                elif message.startswith("Stabilization round "):
                    stab_round = int(message.split()[-1])

                elif message.startswith("Pricing prob "):
                    pricing_prob = int(message.split()[2])

                    # check if the pricing prob should be included in the plot
                    if params['root_only'] and node > 1:
                        continue

                    # store all indices
                    ind_node.append(node)
                    ind_pricing_round.append(pricing_round)
                    ind_stab_round.append(stab_round)
                    ind_pricing_prob.append(pricing_prob)

                    # store the data
                    val_time.append(float(message.split()[-1]))
                    val_nVars.append(int(message.split()[5]))

def main():
    """Entry point when calling this script"""
    args = sys.argv[1:]
    parsed_args = parse_arguments(args)
    set_params(parsed_args)
    if not os.path.exists(params['outdir']):
        os.makedirs(params['outdir'])
    generate_files(parsed_args.filename)

# Calling main script
if __name__ == '__main__':
    main()
