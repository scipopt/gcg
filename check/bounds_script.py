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

def generate_files(files):
    """
    Parse the files and generate temporary files containing only problem name and execution time
    :param files: List of files to be parsed
    :return: A list of all the generated files to be deleted after performance profiling
    """
    for file in files:
        # file = os.path.join(DIR, filename)
        with open(file) as _file:
            df=pd.DataFrame()
            orig = False
            name = None
            read = False
            settings = 'default'
            for line in _file:
                if line.startswith("loaded parameter file"):
                    settings=line.split()[-1]
                    #settings=settings[3]
                    settings=settings.split("/")[-1]
                    settings = os.path.splitext(settings)[0]
                elif line.startswith("Original Program statistics:"):
                    orig = True
                elif line.startswith("Master Program statistics:"):
                    orig = False
                elif line.startswith("Presolved Problem  :"):
                    orig = False
                elif orig and line.startswith("  Problem name     :"):
                    name = line.split()[3]      
                    name = name.split("/")[-1]
                    tmp_name = name.split(".")[-1]
                    if tmp_name[-1] == "gz" or tmp_name[-1] == "z" or tmp_name[-1] == "GZ" or tmp_name[-1] == "Z":
                        name = os.path.splitext(name)[0]
                    name = os.path.splitext(name)[0]
                    print name
                elif line.startswith("Root bounds"):
                    read = True
                elif line.startswith("iter	pb	db"):
                    line_array = line.split()
                    df = pd.DataFrame(columns = line_array, dtype = float)
                elif line.startswith("Pricing Summary:"):
                    read = False
                    df = df.set_index('iter')

                    print df
                    if df.empty:
                        continue
                    
                    df=df.astype(float)
                    
                    #fig, axes = plt.subplots(nrows=2, ncols=2)
                    ax = None
                    ax = df.plot(kind='bar', y='dualdiff', color='green', label='dualdiff', ax=ax, secondary_y=True, alpha=0.5);
                    ax = plt.gca()
                    myLocator = mticker.MultipleLocator(10 ** (math.floor(math.log10(len(df.index)))))
                    ax.xaxis.set_major_locator(myLocator)
                    ax = df.plot(kind='line', y='pb', color='red', label='pb', ax=ax);
                    ax = df.plot(kind='line', y='db', color='blue', label='db', ax=ax);
                    

                    #fig = plot.get_figure()
                    plt.savefig(params['outdir']+"/"+name+"_"+settings+".png")
                elif read:
                    line_array = line.split()
                    df.loc[line_array[0]] = line_array
                               

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
