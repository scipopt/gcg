#!/usr/bin/env python2

import sys
import os
import argparse
from subprocess import call

import glob
import re
import subprocess

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import cm

import matplotlib.ticker as mticker

import math

from matplotlib import gridspec





DIR, FILENAME = os.path.split(__file__)
params = {}
TEMPNAME = 'temp_{}'


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

    parser.add_argument('-x', '--xaxis', type=str,
                        default="time",
                        help='Values to be used in x-axis (can be "time" or "iter"; default: "time")')

    parser.add_argument('filename', nargs='+',
                        help='Name of the files to be used for the creating the bound plots')
    parsed_args = parser.parse_args(args)
    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['outdir'] = args.outdir
    params['xaxis'] = args.xaxis

def generate_files(files):
    """
    Parse the files and generate temporary files containing only problem name and execution time
    :param files: List of files to be parsed
    :return: A list of all the generated files to be deleted after performance profiling
    """
    xaxis = params['xaxis']
    for file in files:
        # file = os.path.join(DIR, filename)
        with open(file) as _file:
            df=pd.DataFrame()
            dfvar=pd.DataFrame()
            orig = False
            name = None
            rootbounds = False
            vardetails = False
            settings = 'default'
            varlines = {}
            varheader = None
            boundlines = {}
            boundheader = None
            for line in _file:
                if line.startswith("loaded parameter file"):
                    # store current settings
                    settings=line.split()[-1]
                    settings=settings.split("/")[-1]
                    settings = os.path.splitext(settings)[0]
                elif not orig and line.startswith("Original Program statistics:"):
                    orig = True
                elif orig and line.startswith("Master Program statistics:"):
                    orig = False
                elif orig and line.startswith("Presolved Problem  :"):
                    orig = False
                elif orig and line.startswith("  Problem name     :"):
                    # store problem name
                    name = line.split()[3]
                    name = name.split("/")[-1]
                    tmp_name = name.split(".")[-1]
                    if tmp_name[-1] == "gz" or tmp_name[-1] == "z" or tmp_name[-1] == "GZ" or tmp_name[-1] == "Z":
                        name = os.path.splitext(name)[0]
                    name = os.path.splitext(name)[0]
                    print name
                elif not rootbounds and line.startswith("Root bounds"):
                    # prepare storage of root bounds
                    rootbounds = True
                elif rootbounds and line.startswith("iter	pb	db"):
                    # store root bounds header
                    line_array = line.split()
                    boundheader = line_array
                elif rootbounds and line.startswith("Pricing Summary:"):
                    # finished with storing root bounds
                    rootbounds = False
                elif rootbounds:
                    # store root bound line
                    line_array = line.split()
                    boundlines[line_array[0]] = line_array
                elif not vardetails and line.startswith("AddedVarDetails:"):
                    # prepare storage of var details
                    vardetails = True
                elif vardetails and line.startswith("VAR: name	node	time") and vardetails:
                    # store var details header
                    line_array = line.split()
                    varheader = line_array[1:]
                elif vardetails and line.startswith("VAR:") and not int(line.split()[2]) == 1:
                    # ignore variables that were not create in the root node
                    continue
                elif vardetails and line.startswith("Root node:"):
                    # finished reading var details (and root bounds)
                    vardetails = False

                    # create dict with root bounds header
                    boundmap = {}
                    for i in range(len(boundheader)):
                        boundmap[i] = boundheader[i]

                    # use boundlines dict to create data frame
                    df = pd.DataFrame.from_dict(data = boundlines, orient = 'index', dtype = float)

                    # if no root bounds are present, ignore instance
                    if len(df) == 0:
                        print "   -> ignored"
                        continue

                    # use root bounds header to rename columns of data frame
                    df.rename(columns = boundmap, inplace=True)

                    # sort lines according to iteration
                    df.sort_values(by='iter', inplace=True)

                    # create var data frame from varlines dict
                    dfvar = pd.DataFrame.from_dict(data = varlines, orient = 'index', dtype = float)

                    # create dict with var header
                    varmap = {}
                    for i in range(len(varheader)):
                        varmap[i] = varheader[i]

                    # use var header to rename columns of var data frame
                    dfvar.rename(columns = varmap, inplace=True)

                    # set index of var data frame to name of var
                    dfvar = dfvar.set_index(keys='name')

                    # set type of var data frame
                    dfvar=dfvar.astype(float)

                    # create new column in data frame containing the number of lp vars generated in each iteration
                    df['nlpvars'] = 0
                    for i in range(len(df)):
                        df.set_value(str(i), 'nlpvars', len(dfvar[(dfvar['rootlpsolval'] > 0) & (dfvar['rootredcostcall'] == i)]))

                    # create new column in data frame containing the number of lp vars generated until each iteration
                    df['nlpvars_cum'] = df[(df['iter'] <= i)].cumsum(axis=0)['nlpvars']

                    # compute total number of vars in root lp solution
                    nlpvars_total = len(dfvar[dfvar['rootlpsolval'] > 0])

                    # create new column in data frame containing the percentage of lp vars generated until each iteration
                    df['lpvars'] = df['nlpvars_cum']/nlpvars_total

                    # repeat this for the vars (generated at the root) in ip solution
                    df['nipvars'] = 0

                    for i in range(len(df)):
                        df.set_value(str(i), 'nipvars', len(dfvar[(dfvar['solval'] > 0) & (dfvar['rootredcostcall'] == i)]))

                    df['nipvars_cum'] = df[(df['iter'] <= i)].cumsum(axis=0)['nipvars']

                    nipvars_total = len(dfvar[dfvar['solval'] > 0])

                    df['ipvars'] = df['nipvars_cum']/nipvars_total

                    # set type of data frame
                    df=df.astype(float)

                    # set infty
                    infty = 10.0 ** 20

                    # set dual bounds of -infinity to NAN
                    df['db'][df['db'] <= -infty] = np.nan

                    # workaround for iterations that were done at the same time (SCIP only counts in 1/100 of a second)
                    df['time_count'] = df.groupby('time')['time'].transform('count')
                    df['time_first'] = df.groupby('time')['iter'].transform('first')
                    df['time'] = df['time'] + 0.01*(df['iter'] - df['time_first'])/df['time_count']
                    df['time_diff'] = df["time"].diff(1)
                    df['time_diff'][0] = df['time'][0]
                    
                    df['db_ma'] = df['db'].rolling(window=5,center=False).mean()

                    # set maximum and minimum of x values (time or iterations) to synchronize the plots
                    xmax = df[xaxis].max()
                    xmin = df[xaxis].min()

                    # set index to time or oterations (depending on which is used)
                    df = df.set_index(keys=xaxis, drop=False)

                    # create grid of 3 plots
                    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
                    ax = plt.subplot(gs[0])
                    ax1 = plt.subplot(gs[1])
                    ax2 = plt.subplot(gs[2])

                    # lp vars plot
                    ax1.set_ylim(bottom=0.0, top=1.1)
                    ax1.set_xlim(left=xmin, right=xmax)
                    ax1 = df.plot(kind='scatter', x=xaxis, y='lpvars', color='blue', label='lpvars', ax=ax1, secondary_y=False, s=1);
                    ax1.set_xticklabels([])
                    x_axis = ax1.axes.get_xaxis()
                    x_axis.set_label_text('')
                    x_axis.set_visible(False)

                    # ip vars plot
                    ax2.set_ylim(bottom=0.0, top=1.1)
                    ax2.set_xlim(left=xmin, right=xmax)
                    ax2 = df.plot(kind='scatter', x=xaxis, y='ipvars', color='red', label='ipvars', ax=ax2, secondary_y=False, s=1);

                    # set base for x labels
                    if( xmax > 0 ):
                        base = 10.0 ** (math.floor(math.log10(xmax)))
                    else:
                        base = 0.01
                    myLocator = mticker.MultipleLocator(base)

                    # specify labels etc. of plot
                    if(xaxis == 'iter' or base > 0.5):
                        majorFormatter = mticker.FormatStrFormatter('%d')
                    else:
                        majorFormatter = mticker.FormatStrFormatter('%0.2f')
                    ax2.xaxis.set_major_locator(myLocator)
                    ax2.xaxis.set_major_formatter(majorFormatter)
                    fixedFormatter = mticker.FormatStrFormatter('%g')
                    ax2.xaxis.set_major_formatter(fixedFormatter)
                    ax2.xaxis.set_minor_locator(plt.NullLocator())
                    lim = ax2.get_xlim()
                    ax2.set_xticks(list(ax2.get_xticks()) + [xmax])
                    ax2.set_xlim(lim)
                    ax2.xaxis.get_major_ticks()[-1].set_pad(15)

                    # set limits and lables for bounds/dualdiff  plot
                    ax.set_xticklabels([])
                    x_axis = ax.axes.get_xaxis()
                    x_axis.set_label_text('')
                    x_axis.set_visible(False)
                    ax.set_xlim(left=xmin, right=xmax)

                    # bounds/dualdiff plot
                    ax = df.plot(kind='line', y='pb', color='red', label='pb', ax=ax, linewidth=0.5);
                    ax = df.plot(kind='line', y='db', color='blue', label='db', ax=ax, linewidth=0.5);
                    ax = df.plot(kind='line', y='db_ma', color='purple', label='db', ax=ax, linewidth=0.5);
                    ax = df.plot(kind='scatter', x=xaxis, y='db', color='blue', label=None, ax=ax, s=0.5);
                    ax = df.plot(kind='line', y='dualdiff', color='green', label='dualdiff', ax=ax, secondary_y=True, alpha=0.25, linewidth=1);
                    ax = df.plot(kind='line', y='dualoptdiff', color='orange', label='dualoptdiff', ax=ax, secondary_y=True, alpha=0.25, linewidth=1);

                    # set y label of secondary y-axis
                    plt.ylabel('diff', fontsize=10, rotation=-90, labelpad=15)

                    # save figure
                    plt.savefig(params['outdir']+"/"+name+"_"+settings+"_"+xaxis+".png")

                    # reset python variables for next instance
                    df = None
                    dfvar = None
                    boundheader = None
                    varheader = None
                    varlines = {}
                    boundlines = {}

                    print "   -> success"
                elif vardetails:
                    # store details of variable
                    line_array = line.split()
                    varlines[line_array[1]] = line_array[1:]

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
