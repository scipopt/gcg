#!/usr/bin/env python3.4

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
                        help='Arguments to be passed on to the performance profiler')
                        
    parser.add_argument('-x', '--xaxis', type=str,
                        default="time",
                        help='Values to be used in x-axis (can be "time" or "iter")')

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
                    #print line
                    name = line.split()[3]      
                    name = name.split("/")[-1]
                    tmp_name = name.split(".")[-1]
                    if tmp_name[-1] == "gz" or tmp_name[-1] == "z" or tmp_name[-1] == "GZ" or tmp_name[-1] == "Z":
                        name = os.path.splitext(name)[0]
                    name = os.path.splitext(name)[0]
                    print name
                elif not rootbounds and line.startswith("Root bounds"):
                    rootbounds = True
                elif rootbounds and line.startswith("iter	pb	db") and rootbounds:
                    line_array = line.split()
                    #df = pd.DataFrame(columns = line_array, dtype = float)
                    boundheader = line_array
                elif rootbounds and line.startswith("Pricing Summary:") and rootbounds:
                    rootbounds = False
                elif rootbounds:
                    line_array = line.split()
                    #df.loc[line_array[0]] = line_array
                    boundlines[line_array[0]] = line_array
                elif not vardetails and line.startswith("AddedVarDetails:"):
                    vardetails = True
                elif vardetails and line.startswith("VAR: name	node	time") and vardetails:
                    line_array = line.split()
                    #dfvar = pd.DataFrame(columns = line_array[1:], dtype = float)
                    varheader = line_array[1:]
                elif vardetails and line.startswith("VAR:") and not int(line.split()[2]) == 1:
                    continue
                elif vardetails and line.startswith("Root node:"):
                    vardetails = False
                                     
                    boundmap = {}
                    for i in range(len(boundheader)):
                        boundmap[i] = boundheader[i]
                       
                    df = pd.DataFrame.from_dict(data = boundlines, orient = 'index', dtype = float)
                    
                    if len(df) == 0:
                        continue
                    #df.colums = boundheader
                    df.rename(columns = boundmap, inplace=True)

                    df.sort_values(by='iter', inplace=True)
                    
                    dfvar = pd.DataFrame.from_dict(data = varlines, orient = 'index', dtype = float)

                    varmap = {}
                    for i in range(len(varheader)):
                        varmap[i] = varheader[i]

                    dfvar.rename(columns = varmap, inplace=True)
                    
                    #print dfvar

                    dfvar = dfvar.set_index(keys='name')

                    dfvar=dfvar.astype(float)

                    df['nlpvars'] = 0

                    for i in range(len(df)):
                        df.set_value(str(i), 'nlpvars', len(dfvar[(dfvar['rootlpsolval'] > 0) & (dfvar['rootredcostcall'] == i)]))

                    df['nipvars'] = 0

                    for i in range(len(df)):
                        df.set_value(str(i), 'nipvars', len(dfvar[(dfvar['solval'] > 0) & (dfvar['rootredcostcall'] == i)]))
                    
                    if df.empty:
                        continue
                    
                    df=df.astype(float)
                    
                    infty = 10.0 ** 20
                    df['db'][df['db'] <= -infty] = np.nan                    
                    
                    lpmax = df['nlpvars'].max()
                    ipmax = df['nipvars'].max()
                    
                    if( np.isnan(ipmax) ):
                        ipmax = 0.0    

                    if( np.isnan(lpmax) ):
                        lpmax = 0.0    
                        
                    #xmax = len(df) - 1
                    
                    #print df['time'].value_counts()
                    #df_size = df.groupby('time').size() 
                    
                    #df_first =  df.groupby('time').first() 

                    df['time_count'] = df.groupby('time')['time'].transform('count') 
                    
                    df['time_first'] = df.groupby('time')['iter'].transform('first') 
                    
                    df['time'] = df['time'] + 0.01*(df['iter'] - df['time_first'])/df['time_count']
                    
                    
                    df['time_diff'] = df["time"].diff(1)
                    df['time_diff'][0] = df['time'][0]

                    xmax = df[xaxis].max()
                    xmin = df[xaxis].min()
                                    
                    df = df.set_index(keys=xaxis, drop=False)
                    
                    #print df
                    
                    #gs = gridspec.GridSpec(4, 1, height_ratios=[3, 1, 1, 1]) 
                    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1]) 
                    ax = plt.subplot(gs[0])
                    ax1 = plt.subplot(gs[1])
                    ax2 = plt.subplot(gs[2])
                    #ax3 = plt.subplot(gs[3])
                    
                    ax1.set_ylim(bottom=0.0, top=lpmax+1) 
                    ax1.set_xlim(left=xmin, right=xmax)
                    ax1 = df.plot(kind='scatter', x=xaxis, y='nlpvars', color='blue', label='nlpvars', ax=ax1, secondary_y=False, s=6);
                    ax1.set_xticklabels([])
                    x_axis = ax1.axes.get_xaxis()
                    x_axis.set_label_text('')
                    x_axis.set_visible(False)                    
                    #base = 10 ** (math.floor(math.log10(len(df.index))))
                    if( xmax > 0 ):
                        base = 10.0 ** (math.floor(math.log10(xmax)))
                    else:
                        base = 0.01
                    ax2.set_ylim(bottom=0.0, top=ipmax+1)
                    ax2.set_xlim(left=xmin, right=xmax)
                    ax2 = df.plot(kind='scatter', x=xaxis, y='nipvars', color='red', label='nipvars', ax=ax2, secondary_y=False, s=6);
                    myLocator = mticker.MultipleLocator(base)
                    if(xaxis == 'iter' or base > 0.5):
                        majorFormatter = mticker.FormatStrFormatter('%d')
                    else:
                        majorFormatter = mticker.FormatStrFormatter('%0.2f')
                    ax2.xaxis.set_major_locator(myLocator)
                    ax2.xaxis.set_major_formatter(majorFormatter)
                    ax2.xaxis.set_minor_locator(plt.NullLocator())                
                    
                    #ax3.set_ylim(bottom=0.0, top=ipmax)
                    #ax3.set_xlim(left=xmin, right=xmax)
                    ##ax3 = df.plot(kind='scatter', x=xaxis, y='time_diff', color='green', label='ptime', ax=ax3, secondary_y=False, s=6);
                    #ax3 = df.plot(kind='line', x=xaxis, y='time_diff', color='green', label='ptime', ax=ax3, secondary_y=False);
                    #myLocator = mticker.MultipleLocator(base)
                    #if(xaxis == 'iter' or base > 0.5):
                    #    majorFormatter = mticker.FormatStrFormatter('%d')
                    #else:
                    #    majorFormatter = mticker.FormatStrFormatter('%0.2f')
                    #ax3.xaxis.set_major_locator(myLocator)
                    #ax3.xaxis.set_major_formatter(majorFormatter)
                    #ax3.xaxis.set_minor_locator(plt.NullLocator())                

                    
                    ax.set_xticklabels([])
                    x_axis = ax.axes.get_xaxis()
                    x_axis.set_label_text('')
                    x_axis.set_visible(False)
                    
                    ax.set_xlim(left=xmin, right=xmax)

                    ax = df.plot(kind='line', y='pb', color='red', label='pb', ax=ax, linewidth=0.5);
                    ax = df.plot(kind='line', y='db', color='blue', label='db', ax=ax, linewidth=0.5);
                    ax = df.plot(kind='scatter', x=xaxis, y='db', color='blue', label=None, ax=ax, s=0.5);
                    ax = df.plot(kind='line', y='dualdiff', color='green', label='dualdiff', ax=ax, secondary_y=True, alpha=0.25, linewidth=1);
                    ax = df.plot(kind='line', y='dualoptdiff', color='orange', label='dualoptdiff', ax=ax, secondary_y=True, alpha=0.25, linewidth=1);
                    
                    plt.savefig(params['outdir']+"/"+name+"_"+settings+"_"+xaxis+".png")
                    
                    df = None
                    dfvar = None
                    
                    boundheader = None
                    varheader = None
                    varlines = {}
                    boundlines = {}
                    
                elif vardetails:
                    line_array = line.split()
                    varlines[line_array[1]] = line_array[1:]
                    #dfvar.loc[line_array[1]] = line_array[1:]


                    
                               

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
