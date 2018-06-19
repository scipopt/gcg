#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 12:25:35 2018

@author: witt
"""
#!/usr/bin/python
import matplotlib  as mpl
mpl.use('Agg')

import os
import sys
import glob
import re
import argparse
import subprocess
#from subprocess import call


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages

DIR, FILENAME = os.path.split(__file__)
params = {}

#os.chdir(os.path.dirname(os.path.realpath(__file__)))

def parse_arguments(args):
    """
    Parse the command-line arguments
    :param args: Command-line arguments provided during execution
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-tm', '--time', type=int, default=20,
                        help='Column for solution time res-files')
    parser.add_argument('-st', '--status', type=int, default=21,
                        help='Column for status in res-file')
    parser.add_argument('-m', '--min', type=float, default=1.0,
                        help='Minimum solution times (solution times are set to maximum of solution times and this value)')
    parser.add_argument('-ss', '--stepsize', type=float, default=0.01,
                        help='Step size used in plot (note that creating the plot can take long when stepsize is too small)')

    parser.add_argument('-l', '--log', type=bool, default=False,
                        help='Should a logarithmic scale be used?')

    parser.add_argument('-o', '--out', type=str, default='perprof_plot.pdf',
                        help='Name of out file')

    parser.add_argument('resfiles', nargs='+',
                        help='Name of .res files to be used for the plot (only instances solved in all .res files are used)')

    parsed_args = parser.parse_args(args)
    return parsed_args


def set_params(args):
    """
    Set the global parame<rs from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['time'] = args.time
    params['status'] = args.status
    params['min'] = args.min
    params['stepsize'] = args.stepsize
    params['log'] = args.log
    params['out'] = args.out

def main():
    """Entry point when calling this script"""
    args = sys.argv[1:]
    parsed_args = parse_arguments(args)
    set_params(parsed_args)

    df=pd.DataFrame()
    resfiles = parsed_args.resfiles

    for resfile in resfiles:
        files = filter(os.path.isfile, glob.glob(resfile))
        files.sort(key=lambda x: os.path.getmtime(x))
        dftmp = pd.DataFrame()
        if len(files) > 0:
            print "res file    :", files[0]

            dftmp = pd.read_csv(
                    files[0],
                    skiprows=3, index_col = 0, delim_whitespace=True, skipfooter=10, usecols=[0,params['time'],params['status']], names=['instance',(resfile,'time'),(resfile,'status')], engine='python')
        else:
            failed = True
            print "WARNING      : No res-file", resfile
        #df = pd.merge(left=dftmp,right=df,how='left',left_index=True, right_index=True)
        df = pd.concat([dftmp,df], axis=1, join='outer')

    for resfile in resfiles:
        df = df[(df[resfile,"status"] == "ok") | (df[resfile,"status"] == "solved")]
        #| (df[resfile,"status"] == "nodelimit")]
        df[resfile,'time'] = df[resfile,'time'].apply(lambda x: max(x, params['min']))

    df["best"] = df[[(resfile, 'time') for resfile in resfiles]].min(axis=1)

    for resfile in resfiles:
        df[resfile,"performance"] = df[resfile,"time"]/df["best"]

    df["rmax"] = df[[(resfile, 'performance') for resfile in resfiles]].max(axis=1)
    rmax = df["rmax"].max()

    def perplot(resfile, t):
        return 1.0*len(df[df[resfile,"performance"] <= t]) / len(df)

    xx = np.arange(1.0,rmax+params['stepsize'],params['stepsize'])

    f, ax = plt.subplots()

    for resfile in resfiles:
        setting = resfile.split(".")[-3]

        yy = np.array([perplot(resfile,x) for x in xx])
        if params['log']:
            ax.semilogx(xx,yy,label=setting)
        else:
            ax.plot(xx,yy,label=setting)

#        ax = df.plot(kind='line',y=(resfile,"per"), ax=ax);
    plt.ylim([0.0, 1.01])
    plt.legend(loc='lower right')

    pp = PdfPages(params['out'])
    plt.savefig(pp, format='pdf')
    pp.close()


# Calling main script
if __name__ == '__main__':
    main()































