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
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-x', '--xaxis', type=str,
                        default="time",
                        choices=['time','iter'],
                        help='Values to be used in x-axis (can be "time" or "iter"; default: "time")')

    parser.add_argument('-f', '--farkas', action='store_true',
                        default=False,
                        help='Include variables produced by Farkas-Pricing in the plot')

    parser.add_argument('-dd', '--dualdiff', action='store_true',
                        default=False,
                        help='Plot difference from current to last dual solution')

    parser.add_argument('-dod', '--dualoptdiff', action='store_true',
                        default=False,
                        help='Plot difference from current to optimal dual solution')

    parser.add_argument('-a', '--average', action='store_true',
                        default=False,
                        help='Plot the moving average for the bounds')

    parser.add_argument('-nobd', '--nobounds', action='store_true',
                        default=False,
                        help='Disable the bounds-subplot')

    parser.add_argument('-bdl', '--boundslinestyle', type=str,
                        default="line",
                        choices=['line','scatter','both'],
                        help='Linestyle of the bounds-plot (can be "line" or "scatter" or "both"; default "line")')

    parser.add_argument('-nolp', '--nolpvars', action='store_true',
                        default=False,
                        help='Disable the lpvars-subplot')

    parser.add_argument('-lpl', '--lplinestyle', type=str,
                        default='line',
                        choices=['line','scatter'],
                        help='Linestyle of the lpvars-plot (can be "line" or "scatter"; default "line")')

    parser.add_argument('-noip', '--noipvars', action='store_true',
                        default=False,
                        help='Disable the ipvars-subplot')

    parser.add_argument('-ipl', '--iplinestyle', type=str,
                        default="line",
                        choices=['line','scatter'],
                        help='Linestyle of the ipvars-plot (can be "line" or "scatter"; default "line")')

    parser.add_argument('-nocmp', '--nocompare', action='store_true',
                        default=False,
                        help='Disable the comparison of different runs of one instance')

    parser.add_argument('filename', nargs='+',
                        help='Names of the files to be used for creating the bound plots')
    parsed_args = parser.parse_args(args)

    # check that at least one subplot is enabled, otherwise exit
    if parsed_args.nobounds and parsed_args.nolpvars and parsed_args.noipvars:
        print 'All plots are disabled. Exiting script.'
        exit()

    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['outdir'] = args.outdir
    params['xaxis'] = args.xaxis
    params['farkas'] = args.farkas
    params['dualdiff'] = args.dualdiff
    params['dualoptdiff'] = args.dualoptdiff
    params['average'] = args.average
    params['bounds'] = not args.nobounds
    params['bdlinestyle'] = args.boundslinestyle
    params['lpvars'] = not args.nolpvars
    params['lplinestyle'] = args.lplinestyle
    params['ipvars'] = not args.noipvars
    params['iplinestyle'] = args.iplinestyle
    params['compare'] = not args.nocompare

def generate_files(files):
    """
    Parse the files and generate temporary files containing only problem name and execution time
    :param files: List of files to be parsed
    :return: A list of all the generated files to be deleted after performance profiling
    """
    xaxis = params['xaxis']

    # Create a dictionary, where all the dataframes, that are generated in the following, are 'globally' stored for comparison
    df_dict = {}
    set_dict = {}

    for file in files:
        # file = os.path.join(DIR, filename)
        with open(file) as _file:
            df = None
            dfvar = None
            orig = False
            name = None
            problemFileName = None
            rootbounds = False
            vardetails = False
            settings = 'default'
            varlines = {}
            varheader = None
            boundlines = {}
            boundheader = None
            scip_status = ""
            for line in _file:
                if line.startswith("@01"):
                    # reset python variables for next instance
                    problemFileName = None
                    df = None
                    dfvar = None
                    boundheader = None
                    varheader = None
                    varlines = {}
                    boundlines = {}
                if line.startswith("loaded parameter file"):
                    # store current settings
                    settings=line.split()[-1]
                    settings=settings.split("/")[-1]
                    settings = os.path.splitext(settings)[0]
                elif not problemFileName and line.startswith("read problem "):
                    # get the problem name from the file name as in "check.awk", in case it is "BLANK" in the actual "Problem name"-line
                    tmparray = line.split("<")[-1].replace(">","").replace("\n","").split("/")[-1].split(".")
                    problemFileName = tmparray[0]
                    if tmparray[-1] == "gz" or tmparray[-1] == "z" or tmparray[-1] == "GZ" or tmparray[-1] == "Z":
                        tmparray.pop()
                    for i in range(1,len(tmparray)-1):
                        problemFileName += "." + tmparray[i]
                elif not orig and line.startswith("Original Program statistics:"):
                    orig = True
                elif orig and line.startswith("Master Program statistics:"):
                    orig = False
                elif orig and line.startswith("Presolved Problem  :"):
                    orig = False
                elif orig and line.startswith("SCIP Status        :"):
                    scip_status = line.split(":")[-1].strip()
                elif orig and line.startswith("  Problem name     :"):
                    # store problem name
                    name = line.split()[3]
                    name = name.split("/")[-1]
                    tmp_name = name.split(".")[-1]
                    if tmp_name[-1] == "gz" or tmp_name[-1] == "z" or tmp_name[-1] == "GZ" or tmp_name[-1] == "Z":
                        name = os.path.splitext(name)[0]
                    name = os.path.splitext(name)[0]
                    if name == 'BLANK':
                        name = problemFileName
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
                        print name
                        print "   -> ignored"
                        print "   -> SCIP Status : " + scip_status
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
                        df.set_value(str(i), 'nlpvars', len(dfvar[(dfvar['rootlpsolval'] <> 0) & (dfvar['rootredcostcall'] == i)]))

                    # add the number of all lp-variables, not created by reduced cost pricing (e.g. by Farkas-Pricing)
                    if params['farkas']:
                        df.set_value(str(0),'nlpvars', df['nlpvars'][0] + len(dfvar[(dfvar['rootlpsolval'] <> 0) & (dfvar['rootredcostcall'] == -1.)]))

                    # create new column in data frame containing the number of lp vars generated until each iteration
                    df['nlpvars_cum'] = df[(df['iter'] < len(df))].cumsum(axis=0)['nlpvars']

                    # compute total number of vars in root lp solution, that are included in the plot
                    nlpvars_total = df['nlpvars_cum'].iloc[-1]

                    # create new column in data frame containing the percentage of lp vars generated until each iteration
                    df['lpvars'] = df['nlpvars_cum']/nlpvars_total

                    # repeat this for the vars (generated at the root) in ip solution
                    df['nipvars'] = 0

                    for i in range(len(df)):
                        df.set_value(str(i), 'nipvars', len(dfvar[(dfvar['solval'] > 0) & (dfvar['rootredcostcall'] == i)]))

                    if params['farkas']:
                        df.set_value(str(0),'nipvars', df['nipvars'][0] + len(dfvar[(dfvar['solval'] > 0) & (dfvar['rootredcostcall'] == -1.)]))

                    df['nipvars_cum'] = df[(df['iter'] < len(df))].cumsum(axis=0)['nipvars']

                    nipvars_total = df['nipvars_cum'].iloc[-1]

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
                    # workaround for identical limits (will produce UserWarnings otherwise)
                    if xmax == xmin:
                        if xmax == 0:
                            xmax = 0.01
                        else:
                            xmin = 0.95*xmin
                            xmax = 1.05*xmax

                    # set index to time or iterations (depending on which is used)
                    df = df.set_index(keys=xaxis, drop=False)

                    # number of plots, the user wants
                    nplots = params['bounds'] + params['lpvars'] + params['ipvars']

                    # create grid of nplots plots
                    if params['bounds']:
                        height_ratios = [3]+[1]*(nplots-1)
                    else:
                        height_ratios = [1]*nplots
                    gs = list(gridspec.GridSpec(nplots, 1, height_ratios=height_ratios))
                    axes = {}
                    if params['ipvars']:
                        axes['ip'] = plt.subplot(gs.pop())
                    if params['lpvars']:
                        axes['lp'] = plt.subplot(gs.pop())
                    if params['bounds']:
                        axes['db'] = plt.subplot(gs.pop())

                    # lp vars plot
                    if params['lpvars']:
                        axes['lp'].set_ylim(bottom=0.0, top=1.1)
                        axes['lp'].set_xlim(left=xmin, right=xmax)
                        frmtStr = 'c'
                        if params['lplinestyle'] == 'line':
                            frmtStr += '-'
                        elif params['lplinestyle'] == 'scatter':
                            frmtStr += 'o'
                        axes['lp'].plot(df[xaxis], df['lpvars'], frmtStr, label ='lpvars', markersize=1.6, linewidth = 0.8)
                        axes['lp'].set_ylabel('lpvars')
                        axes['lp'].set_xticklabels([])
                        x_axis = axes['lp'].axes.get_xaxis()
                        x_axis.set_label_text('')
                        x_axis.set_visible(False)

                    # ip vars plot
                    if params['ipvars']:
                        axes['ip'].set_ylim(bottom=0.0, top=1.1)
                        axes['ip'].set_xlim(left=xmin, right=xmax)
                        frmtStr = 'y'
                        if params['iplinestyle'] == 'line':
                            frmtStr += '-'
                        elif params['iplinestyle'] == 'scatter':
                            frmtStr += 'o'
                        axes['ip'].plot(df[xaxis], df['ipvars'], frmtStr, label ='ipvars', markersize=1.6, linewidth = 0.8)
                        axes['ip'].set_ylabel('ipvars')
                        axes['ip'].set_xticklabels([])
                        x_axis = axes['ip'].axes.get_xaxis()
                        x_axis.set_label_text('')
                        x_axis.set_visible(False)

                    if params['bounds']:
                        # set limits and lables for bounds/dualdiff  plot
                        axes['db'].set_xticklabels([])
                        x_axis = axes['db'].axes.get_xaxis()
                        x_axis.set_label_text('')
                        x_axis.set_visible(False)
                        axes['db'].set_xlim(left=xmin, right=xmax)

                        # create a new axis for the difference-plots, since they need a different y-label
                        if params['dualdiff'] or params['dualoptdiff']:
                            axes['db_diff'] = axes['db'].twinx()

                        # bounds/dualdiff plot
                        if params['bdlinestyle'] == 'line':
                            frmtStr = '-'
                        elif params['bdlinestyle'] == 'scatter':
                            frmtStr = 'o'
                        elif params['bdlinestyle'] == 'both':
                            frmtStr = '-o'
                        axes['db'].plot(df[xaxis], df['pb'], frmtStr, color = 'red', label='pb', linewidth=0.8, markersize = 1.6)
                        axes['db'].plot(df[xaxis], df['db'], frmtStr, color = 'blue', label='db', linewidth=0.8, markersize = 1.6)
                        if params['average']:
                            axes['db'].plot(df[xaxis], df['db_ma'], '-', color = 'purple', label='db (average)', linewidth=0.5)
                        if params['dualdiff']:
                            axes['db_diff'].plot(df[xaxis], df['dualdiff'], 'g-', label='dualdiff', alpha = .25, linewidth=1)
                        if params['dualoptdiff']:
                           axes['db_diff'].plot(df[xaxis], df['dualoptdiff'], '-', color = 'orange', label='dualoptdiff', alpha = .25, linewidth=1)

                        # create the legend and set the primary y-label
                        lines, labels = axes['db'].get_legend_handles_labels()
                        if params['dualdiff'] or params['dualoptdiff']:
                            lines += axes['db_diff'].get_legend_handles_labels()[0]
                            labels += axes['db_diff'].get_legend_handles_labels()[1]
                        axes['db'].legend(lines, labels)
                        axes['db'].set_ylabel('Bounds')

                    # set base for x labels
                    if( xmax > 0 ):
                        base = 10.0 ** (math.floor(math.log10(xmax)))
                    else:
                        base = 0.01
                    myLocator = mticker.MultipleLocator(base)

                    # specify labels etc. of plot
                    if params['ipvars']:
                        lowest_ax = axes['ip']
                    elif params['lpvars']:
                        lowest_ax = axes['lp']
                    elif params['bounds']:
                        lowest_ax = axes['db']
                    if(xaxis == 'iter' or base > 0.5):
                        majorFormatter = mticker.FormatStrFormatter('%d')
                    else:
                        majorFormatter = mticker.FormatStrFormatter('%0.2f')
                    lowest_ax.xaxis.set_major_locator(myLocator)
                    lowest_ax.xaxis.set_major_formatter(majorFormatter)
                    fixedFormatter = mticker.FormatStrFormatter('%g')
                    lowest_ax.xaxis.set_major_formatter(fixedFormatter)
                    lowest_ax.xaxis.set_minor_locator(plt.NullLocator())
                    lim = lowest_ax.get_xlim()
                    xmax_rounded = round(xmax, int(-math.log10(base)))
                    if (xmax_rounded in list(lowest_ax.get_xticks())):
                        xticks = list(lowest_ax.get_xticks())
                    else:
                        xticks = list(lowest_ax.get_xticks()) + [xmax_rounded]
                    lowest_ax.set_xticks(xticks)
                    lowest_ax.set_xlim(lim)
                    lowest_ax.xaxis.get_major_ticks()[-1].set_pad(15)
                    lowest_ax.set_xlabel(xaxis)
                    lowest_ax.xaxis.set_visible(True)

                    # set y label of secondary y-axis if necessary
                    if params['dualdiff'] or params['dualoptdiff']:
                        plt.ylabel('Differences', fontsize=10, rotation=-90, labelpad=15)

                    # ensure, that there is enough space for labels
                    plt.tight_layout()

                    # set the size of the figure (a too small size will lead to too large legends)
                    plt.gcf().set_size_inches(9.33,7)

                    # save figure and ensure, that there are not two files with the same name
                    fig_filename = params['outdir']+"/"+name+"_"+settings+"_"+xaxis
                    i = ""
                    while os.path.isfile(fig_filename + i + ".png"):
                        if i == "":
                            i = "2"
                        else:
                            i = str(int(i)+1)
                    plt.savefig(fig_filename + i + ".png", dpi=300)

                    # store the current dataframe globally
                    if params['compare']:
                        if not (name in df_dict):
                            df_dict[name] = []
                            set_dict[name] = []
                        df_dict[name].append(df.copy())
                        set_dict[name].append(settings)

                    # reset python variables for next instance
                    problemFileName = None
                    df = None
                    dfvar = None
                    boundheader = None
                    varheader = None
                    varlines = {}
                    boundlines = {}

                elif vardetails:
                    # store details of variable
                    line_array = line.split()
                    varlines[line_array[1]] = line_array[1:]

    # compare different runs of one instance
    if params['compare']:
        for name, runs in df_dict.iteritems():
            if len(runs) > 1:

                # set maximum and minimum of x values (time or iterations) to synchronize the plots
                infty = 10.0 ** 20
                xmax = -infty
                xmin = infty
                for run in runs:
                    if run[xaxis].max() > xmax:
                        xmax = run[xaxis].max()
                    if run[xaxis].min() < xmin:
                        xmin = run[xaxis].min()

                # workaround for identical limits
                if xmin == xmax:
                    if xmax == 0:
                        xmax = 0.01
                    else:
                        xmin = 0.95*xmin
                        xmax = 1.05*xmax

                # number of plots, the user wants
                nplots = params['bounds'] + params['lpvars'] + params['ipvars']

                # create grid of nplots plots
                if params['bounds']:
                    height_ratios = [3]+[1]*(nplots-1)
                else:
                    height_ratios = [1]*nplots
                gs = list(gridspec.GridSpec(nplots, 1, height_ratios=height_ratios))
                axes = {}
                if params['ipvars']:
                    axes['ip'] = plt.subplot(gs.pop())
                if params['lpvars']:
                    axes['lp'] = plt.subplot(gs.pop())
                if params['bounds']:
                    axes['db'] = plt.subplot(gs.pop())

                # lp vars plot
                if params['lpvars']:
                    axes['lp'].set_ylim(bottom=0.0, top=1.1)
                    axes['lp'].set_xlim(left=xmin, right=xmax)
                    axes['lp'].set_ylabel('lpvars')
                    axes['lp'].set_xticklabels([])
                    x_axis = axes['lp'].axes.get_xaxis()
                    x_axis.set_label_text('')
                    x_axis.set_visible(False)

                # ip vars plot
                if params['ipvars']:
                    axes['ip'].set_ylim(bottom=0.0, top=1.1)
                    axes['ip'].set_xlim(left=xmin, right=xmax)
                    axes['ip'].set_ylabel('ipvars')
                    axes['ip'].set_xticklabels([])
                    x_axis = axes['ip'].axes.get_xaxis()
                    x_axis.set_label_text('')
                    x_axis.set_visible(False)

                # bounds plot
                if params['bounds']:
                    # set limits and lables for bounds/dualdiff  plot
                    axes['db'].set_xticklabels([])
                    x_axis = axes['db'].axes.get_xaxis()
                    x_axis.set_label_text('')
                    x_axis.set_visible(False)
                    axes['db'].set_xlim(left=xmin, right=xmax)
                    axes['db'].set_ylabel('Bounds')

                # set base for x labels
                if( xmax > 0 ):
                    base = 10.0 ** (math.floor(math.log10(xmax)))
                else:
                    base = 0.01
                myLocator = mticker.MultipleLocator(base)

                # specify labels etc. of plot
                if params['ipvars']:
                    lowest_ax = axes['ip']
                elif params['lpvars']:
                    lowest_ax = axes['lp']
                elif params['bounds']:
                    lowest_ax = axes['db']
                if(xaxis == 'iter' or base > 0.5):
                    majorFormatter = mticker.FormatStrFormatter('%d')
                else:
                    majorFormatter = mticker.FormatStrFormatter('%0.2f')
                lowest_ax.xaxis.set_major_locator(myLocator)
                lowest_ax.xaxis.set_major_formatter(majorFormatter)
                fixedFormatter = mticker.FormatStrFormatter('%g')
                lowest_ax.xaxis.set_major_formatter(fixedFormatter)
                lowest_ax.xaxis.set_minor_locator(plt.NullLocator())
                lim = lowest_ax.get_xlim()
                xmax_rounded = round(xmax, int(-math.log10(base)))
                if (xmax_rounded in list(lowest_ax.get_xticks())):
                    xticks = list(lowest_ax.get_xticks())
                else:
                    xticks = list(lowest_ax.get_xticks()) + [xmax_rounded]
                lowest_ax.set_xticks(xticks)
                lowest_ax.set_xlim(lim)
                lowest_ax.xaxis.get_major_ticks()[-1].set_pad(15)
                lowest_ax.set_xlabel(xaxis)
                lowest_ax.xaxis.set_visible(True)

                # build colormaps for the diagrams
                cmap = {}
                for p in axes:
                    if p == 'db':
                        cmap[p] = plt.cm.get_cmap('jet', max(len(runs), 5))
                    else:
                        cmap[p] = plt.cm.get_cmap('jet', max(len(runs), 5))

                # plot all the runs
                # first, create a list, to store the plot-handles, that have to be inlcuded in the legend
                handles = []
                for iter_run, df in enumerate(runs):
                    # plot the lpvars
                    if params['lpvars']:
                        if params['lplinestyle'] == 'line':
                            frmtStr = '-'
                        elif params['lplinestyle'] == 'scatter':
                            frmtStr = 'o'
                        axes['lp'].plot(df[xaxis], df['lpvars'], frmtStr, color = cmap['lp'](iter_run), label ='lpvars ' + set_dict[name][iter_run], markersize=1.6, linewidth = 0.8)

                    # plot the ipvars
                    if params['ipvars']:
                        if params['iplinestyle'] == 'line':
                            frmtStr = '-'
                        elif params['iplinestyle'] == 'scatter':
                            frmtStr = 'o'
                        axes['ip'].plot(df[xaxis], df['ipvars'], frmtStr, color = cmap['ip'](iter_run), label ='ipvars '+ set_dict[name][iter_run], markersize=1.6, linewidth = 0.8)

                    # bounds/dualdiff plot
                    if params['bounds']:
                        if params['bdlinestyle'] == 'line':
                            frmtStr = '-'
                        elif params['bdlinestyle'] == 'scatter':
                            frmtStr = 'o'
                        elif params['bdlinestyle'] == 'both':
                            frmtStr = '-o'
                        tmp, = axes['db'].plot(df[xaxis], df['pb'], frmtStr, color = cmap['db'](iter_run), label=set_dict[name][iter_run], linewidth=0.8, markersize = 1.6)
                        handles.append(tmp)
                        axes['db'].plot(df[xaxis], df['db'], frmtStr, color = cmap['db'](iter_run), linewidth=0.8, markersize = 1.6)

                if params['dualdiff'] or params['dualoptdiff']:
                    # plot the differences
                    axes['db_diff'] = axes['db'].twinx()
                    for iter_run, df in enumerate(runs):
                        if params['dualdiff']:
                            axes['db_diff'].plot(df[xaxis], df['dualdiff'], '--', color = cmap['db'](iter_run), label='dualdiff ' + set_dict[name][iter_run], linewidth=0.8, markersize = 1.6, alpha = 0.6)
                        if params['dualoptdiff']:
                            axes['db_diff'].plot(df[xaxis], df['dualoptdiff'], '--', color = cmap['db'](iter_run), label='dualoptdiff ' + set_dict[name][iter_run], linewidth=0.8, markersize = 1.6, alpha = 0.6)
                    # set y label of secondary y-axis
                    plt.ylabel('Differences', fontsize=10, rotation=-90, labelpad=15)

                # ensure, that there is enough space for labels
                plt.tight_layout()

                # create the legend
                axes['db'].legend(handles = handles, loc='lower left', bbox_to_anchor = (0,1.02,1,0.2), ncol=3, mode='expand')

                # make room for the legend
                plt.subplots_adjust(top=0.9)

                # set the size of the figure (a too small size will lead to too large legends)
                plt.gcf().set_size_inches(9.33,7)

                # save figure and ensure, that there are not two files with the same name
                fig_filename = params['outdir']+"/"+ "compareRuns_" + name+"_"+xaxis
                i = ""
                while os.path.isfile(fig_filename + i + ".png"):
                    if i == "":
                        i = "1"
                    else:
                        i = str(int(i)+1)
                plt.savefig(fig_filename + i + ".png", dpi = 300)


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
