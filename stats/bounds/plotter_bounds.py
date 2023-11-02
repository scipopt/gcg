#!/usr/bin/env python3

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
from matplotlib.ticker import FuncFormatter

import math

from matplotlib import gridspec
import pickle as pickler

if os.path.isdir("bounds"):
    try: import bounds.parser_bounds
    except: import parser_bounds
else:
    import parser_bounds


DIR, FILENAME = os.path.split(__file__)
params = {}
TEMPNAME = 'temp_{}'
xaxis=""

def parse_arguments(args):
    """
    Parse the command-line arguments
    :param args: Command-line arguments provided during execution
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', type=str,
                        default="./plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-x', '--xaxis', type=str,
                        default="time",
                        choices=['time','iter'],
                        help='Values to be used in x-axis (can be "time" or "iter"; default: "time")')

    parser.add_argument('-f', '--farkas', action='store_true',
                        default=False,
                        help='Include variables produced by Farkas-Pricing in the plot')

    parser.add_argument('-png', action='store_true',
                        default=False,
                        help='Save all non-comparison plots (that do not yet exist) as png.')

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

    parser.add_argument('-nogap', '--nogapplot', action='store_true',
                        default=False,
                        help='Disable the gap-subplot')

    parser.add_argument('-gpl', '--gaplinestyle', type=str,
                        default="x",
                        choices=['x','line','scatter'],
                        help='Linestyle of the gap-plot (can be "line" or "scatter"; default "x")')

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

    parser.add_argument('-load', '--loadpickle', action='store_true',
                        default=False,
                        help='Load a pickle, do not parse outfile. Give a folder containing only the boundsset and boundsdict.')

    parser.add_argument('-save', '--savepickle', action='store_true',
                        default=False,
                        help='Save a pickle, do not generate visualizations. Will be saved into dataframedir.')

    parser.add_argument('filename', nargs='+',
                        help='Names of the files to be used for creating the bound plots')

    parsed_args = parser.parse_args(args)

    # check that at least one subplot is enabled, otherwise exit
    if parsed_args.nobounds and parsed_args.nolpvars and parsed_args.noipvars and parsed_args.nogap:
        print('All plots are disabled. Exiting script.')
        exit()

    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['filename'] = args.filename[0]
    params['outdir'] = args.outdir
    params['xaxis'] = args.xaxis
    params['farkas'] = args.farkas
    params['dualdiff'] = args.dualdiff
    params['dualoptdiff'] = args.dualoptdiff
    params['average'] = args.average
    params['bounds'] = not args.nobounds
    params['bdlinestyle'] = args.boundslinestyle
    params['gap'] = not args.nogapplot
    params['gaplinestyle'] = args.gaplinestyle
    params['lpvars'] = not args.nolpvars
    params['lplinestyle'] = args.lplinestyle
    params['ipvars'] = not args.noipvars
    params['iplinestyle'] = args.iplinestyle
    params['compare'] = not args.nocompare
    params['load'] = args.loadpickle
    params['save'] = args.savepickle
    params['png'] = args.png
    return params

def generate_visu(dir, df_dict = {}, set_dict = {}, params = {}):
    xaxis = params['xaxis']
    files = []
    #df_dict = {}
    #set_dict = {}
    boundsFound = False
    dictFound = False
    setFound = False
    vbcinfoFound = False
    vbcFound = False
    if params['load']:
        try:
            if os.listdir(dir) == []:
                print("Warning: dataframe directory empty.\nTerminating.")
                exit()
        except NotADirectoryError:
            print("Warning: not a directory. Please give a directory (e.g. plots/), where the *.boundsdict.pkl and *.boundsset.pkl are in.\nTerminating.")
            exit()
        for file in os.listdir(dir):
            # Load dataframe for finished instance - maybe a feature for the future, if one wants to have pickles of single instances
            #if file.endswith("bounds.pkl"):
            #    files.append(os.path.join(dir, file))
            #    boundsFound = True
            #    #print(file)
            if file.endswith("boundsdict.pkl"):
                with open(os.path.join(dir, file), 'rb') as handle:
                    df_dict = pickler.load(handle)
                    dictFound = True
            if file.endswith("boundsset.pkl"):
                with open(os.path.join(dir, file), 'rb') as handle:
                    set_dict = pickler.load(handle)
                    setFound = True

    #if params['load'] and not boundsFound:
    #    print("Fatal: *.bounds.pkl not found.\nTerminating.")
    #    exit()
    if params['load'] and not dictFound:
        print("Fatal: *.boundsdict.pkl not found.\nTerminating.")
        exit()
    if params['load'] and not setFound:
        print("Fatal: *.boundsset.pkl not found.\nTerminating.")
        exit()

    for instance in df_dict:
        df = df_dict[instance]
        # if params['load']:
        #    df = pd.read_pickle(params["dataframedir"] + "/" + instance + ".bounds.pkl")
        #else:
        df = df_dict[instance][0]

        # append df to df_dict at that place
        name = str(instance)
        settings = set_dict[name][0]

        # set maximum and minimum of x values (time or iterations) to synchronize the plots
        try:
            xaxis = params['xaxis']
            xmax = df[xaxis].max()
            xmin = df[xaxis].min()
        except KeyError:
            print("   -> skipping {}".format(name))
            continue
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
        nplots = params['bounds'] + params['lpvars'] + params['ipvars'] + params['gap']

        # create grid of nplots plots
        if params['bounds']:
            height_ratios = [1]+[0.3]*(nplots-1)
            if params['gap']:
                height_ratios = [2]+[1]+[0.5]*(nplots-2)
        elif params['gap']:
            height_ratios = [5]+[1]*(nplots-1)
        else:
            height_ratios = [1]*nplots
        gs = list(gridspec.GridSpec(nplots, 1, height_ratios=height_ratios))
        axes = {}
        if params['ipvars']:
            axes['ip'] = plt.subplot(gs.pop())
        if params['lpvars']:
            axes['lp'] = plt.subplot(gs.pop())
        if params['gap']:
            axes['gap'] = plt.subplot(gs.pop())
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

        # ip vars plot
        if params['gap']:
            #axes['gap'].set_ylim(bottom=0.0, top=1.1)
            axes['gap'].set_xlim(left=xmin, right=xmax)
            frmtStr = 'x'
            if params['gaplinestyle'] == 'line':
                frmtStr = '-'
            elif params['gaplinestyle'] == 'scatter':
                frmtStr = 'o'
            axes['gap'].yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
            axes['gap'].plot(df[xaxis], df['gap'], frmtStr, label ='gap', markersize=3.6, linewidth = 0.8)
            axes['gap'].set_ylabel('gap')
            axes['gap'].set_xticklabels([])
            x_axis = axes['gap'].axes.get_xaxis()
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

            # Fill all NaN values in the dual bound with the previous entry
            df['db'].fillna(method='ffill',inplace=True)
            axes['db'].plot(df[xaxis], df['pb'], frmtStr, color = 'red', label='primal bound', linewidth=0.8, markersize = 1.6)
            axes['db'].plot(df[xaxis], df['db'], frmtStr, color = 'blue', label='dual bound', linewidth=0.8, markersize = 1.6)
            if params['average']:
                axes['db'].plot(df[xaxis], df['db_ma'], '-', color = 'purple', label='dual bound (average)', linewidth=0.5)
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
        if params['bounds']:
            lowest_ax = axes['db']
        if params['gap']:
            lowest_ax = axes['gap']
        if params['lpvars']:
            lowest_ax = axes['lp']
        if params['ipvars']:
            lowest_ax = axes['ip']
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
        #plt.tight_layout()

        # set the size of the figure (a too small size will lead to too large legends)
        plt.gcf().set_size_inches(9.33,7)

        # save figure and ensure, that there are not two files with the same name
        fig_filename = params['outdir']+"/"+name+"."+settings+".bounds."+xaxis
        i = ""
        while os.path.isfile(fig_filename + i + ".pdf"):
            if i == "":
                i = "2"
            else:
                i = str(int(i)+1)
        plt.figtext(.5,.93,"Instance: {}".format(name.split('/')[-1]),ha="center",size="14")
        plt.title("Primal/Dual Bound Development in the Root Node")
        if params['interactive']:
            yield instance, plt.gcf()
        elif params['png']:
            plt.savefig(fig_filename + i + ".png", dpi=300)
        else:
            plt.savefig(fig_filename + i + ".pdf", dpi=300)
        plt.close()


        # compare different runs of one instance
    if params['compare']:
        for name, runs in df_dict.items():
            if len(runs) > 1:
                abortRun = False
                # set maximum and minimum of x values (time or iterations) to synchronize the plots
                infty = 10.0 ** 20
                xmax = -infty
                xmin = infty
                for run in runs:
                    try:
                        if run[xaxis].max() > xmax:
                            xmax = run[xaxis].max()
                        if run[xaxis].min() < xmin:
                            xmin = run[xaxis].min()
                    except KeyError:
                        print("Information: Could not synchronize xaxis ({}) in comparison plot for instance {}.".format(xaxis,name))
                        abortRun = True
                        break
                if abortRun:
                    break
                # workaround for identical limits
                if xmin == xmax:
                    if xmax == 0:
                        xmax = 0.01
                    else:
                        xmin = 0.95*xmin
                        xmax = 1.05*xmax

                # number of plots, the user wants
                nplots = params['bounds'] + params['lpvars'] + params['ipvars'] + params['gap']

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
                if params['gap']:
                    axes['gap'] = plt.subplot(gs.pop())
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

                # gap plot
                if params['gap']:
                    #axes['gap'].set_ylim(bottom=0.0, top=1.1)
                    axes['gap'].yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
                    axes['gap'].set_xlim(left=xmin, right=xmax)
                    axes['gap'].set_ylabel('gap')
                    axes['gap'].set_xticklabels([])
                    x_axis = axes['gap'].axes.get_xaxis()
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
                if params['bounds']:
                    lowest_ax = axes['db']
                if params['gap']:
                    lowest_ax = axes['gap']
                if params['lpvars']:
                    lowest_ax = axes['lp']
                if params['ipvars']:
                    lowest_ax = axes['ip']

                if params['ipvars']:
                    highest_ax = axes['ip']
                if params['lpvars']:
                    highest_ax = axes['lp']
                if params['gap']:
                    highest_ax = axes['gap']
                if params['bounds']:
                    highest_ax = axes['db']

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
                try:
                    if (xmax_rounded in list(lowest_ax.get_xticks())):
                        xticks = list(lowest_ax.get_xticks())
                    else:
                        xticks = list(lowest_ax.get_xticks()) + [xmax_rounded]
                except ValueError:
                    print("Information: Could not add some ticks.")
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

                    # plot the ipvars
                    if params['gap']:
                        if params['gaplinestyle'] == 'line':
                            frmtStr = '-'
                        elif params['gaplinestyle'] == 'scatter':
                            frmtStr = 'o'
                        elif params['gaplinestyle'] == 'x':
                            frmtStr = '.'
                        axes['gap'].plot(df[xaxis], df['gap'], frmtStr, color = cmap['gap'](iter_run), label ='gap '+ set_dict[name][iter_run], markersize=0.8, linewidth = 0.8)

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
                        df['db'].fillna(method='ffill',inplace=True)
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
                highest_ax.legend(handles = handles, loc='lower left', bbox_to_anchor = (0,1.02,1,0.2), ncol=4, mode='expand')

                # make room for the legend
                plt.subplots_adjust(top=0.85)

                # set the size of the figure (a too small size will lead to too large legends)
                plt.gcf().set_size_inches(9.33,7)

                # save figure and ensure, that there are not two files with the same name
                fig_filename = params['outdir']+"/"+ name + ".compare"+".bounds_" +xaxis
                i = ""
                while os.path.isfile(fig_filename + i + ".pdf"):
                    if i == "":
                        i = "1"
                    else:
                        i = str(int(i)+1)
                plt.figtext(.5,.95,"Instance: {}".format(name.split('/')[-1]),ha="center",size="14")
                highest_ax.set_title("Comparison of Primal/Dual Bound Development in the Root Node", y=1.2)
                if params['interactive']:
                    yield instance, plt.gcf()
                plt.savefig(fig_filename + i + ".pdf", dpi = 300)
                plt.close()


def main():
    """Entry point when calling this script"""
    args = sys.argv[1:]
    parsed_args = parse_arguments(args)
    set_params(parsed_args)
    if params['save'] and params['load']:
        params['save'] = False
        params['load'] = False

    if not params['load']:
        for file in parsed_args.filename:
            if not os.path.isfile(file):
                print("Fatal: File '{}' could not be found.\nTerminating.".format(file))
                exit()
    df_dict = {}
    set_dict = {}
    if parsed_args.filename[0].endswith(".out"):
        print("Parsing files...")
        df_dict, set_dict = parser_bounds.generate_files(parsed_args.filename,params)
    print("Generating visualizations...")
    if not os.path.exists(params['outdir']):
        os.makedirs(params['outdir'])
    generate_visu(params['filename'], df_dict = df_dict, set_dict = set_dict, params = params)

# Calling main script
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('\nKeyboardInterrupt.\nTerminating.')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
