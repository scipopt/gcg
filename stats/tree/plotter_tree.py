#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import argparse
sys.path.append("./")
import misc.vbc_reader as vbcr

params = {}

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

    parser.add_argument('input', nargs=1,
                        default="",
                        help='input directory or file (default: "dataframes")')

    parser.add_argument('-B', '--bar',
                        action="store_true",
                        help='create barchart')

    parser.add_argument('-P', '--plot',
                        action="store_true",
                        help='create simple plot')

    parser.add_argument('-a', '--absolute',
                        action="store_true",
                        help='deactivate normalization of the data (for each depth level, divide through all nodes on this level)')

    parser.add_argument('-save', '--savepickle', action='store_true',
                        help='saves the collected data in a pickle-file, in the OUTDIR, instead of plotting it (see also --loadpickle)')

    parser.add_argument('-load', '--loadpickle', action='store_true',
                        help='loads the collected data from a pickle-file (see also --savepickle)')

    parsed_args = parser.parse_args(args)

    return parsed_args

def set_params(args):
    """
    Set the global parameters from the parsed command-line arguments
    :param args: parsed command-line arguments
    :return:
    """
    params['outdir'] = args.outdir
    params['input'] = args.input
    params['save'] = args.savepickle
    params['load'] = args.loadpickle
    if args.bar and args.plot:
        params['type'] = "both"
    elif args.bar:
        params['type'] = "bar"
    elif args.plot:
        params['type'] = 'plot'
    else:
        print("You did not select which plot to generate. Making both.")
        params['type'] = "both"
    params['normalize'] = not args.absolute
    params['interactive'] = False

#def plotSupervisor(dict,name,settings,type="both"):


def plot(dict,name,settings,params,interactive=False):
    if params["type"] == "both":
        p = params.copy()
        p["type"] = "bar"
        plot(dict,name,settings,p)
        p["type"] = "plot"
        plot(dict,name,settings,p)
        return

    #print(settings)

    # create plotting data and plot itself
    plotdata = [ v for v in dict.values() ]
    fig, ax = plt.subplots()

    # make settings for plot
    #ax = plt.figure(figsize=(25, 10))
    ax.tick_params(axis ='both', which ='both', length = 1)
    #if max(plotdata) <= 30:
    #    ax.yaxis.set_ticks(np.arange(0,max(plotdata),1))
    #elif max(plotdata) > 30 and max(plotdata) < 50:
    #    ax.yaxis.set_ticks(np.arange(0,max(plotdata),2))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    nnodes = sum(plotdata)

    if params["type"] == "plot":
        ax.set(xlabel="Tree depth", ylabel="Opened nodes (absolute)")
        ax.plot(plotdata, 'x', markersize=8)
    elif params["type"] == "bar":
        if params['normalize']:
            ax.set(xlabel="Tree depth", ylabel="Opened nodes (ratio)")
            for i in range(len(plotdata)):
                plotdata[i] = plotdata[i] / 2**i
            ax.bar(range(len(plotdata)),plotdata)
        else:
            ax.set(xlabel="Tree depth", ylabel="Opened nodes (absolute)")
            ax.bar(range(len(plotdata)),plotdata)
    else:
        print("No plot selected.")

    # add texts
    plt.figtext(.5,.93,"Instance: {}".format(name.split('/')[-1]),ha="center",size="14")
    plt.figtext(.01,.01,'The total number of opened nodes was ' + str(nnodes) + '.',size='12')
    plt.title("Number of opened nodes on each tree depth level",size="18",va="bottom")

    # save and close
    
    fig.set_size_inches(18.5, 10.5)
    if params['outdir'] != None: 
        print("Saving plot as {}.{}.tree.{}.pdf".format(os.path.join(params['outdir'],name.split('/')[-1]),settings,params["type"]))
        fig.savefig("{}.{}.tree.{}.pdf".format(os.path.join(params['outdir'],name.split('/')[-1]),settings,params["type"]),dpi=300)
    if params['interactive']: 
        return fig

def calc_itperdepth(vbc_df):
    print("Calculating depths")
    it_per_depth = {}
    # fill dict
    for num in vbc_df['depth']:
        it_per_depth[num] = 0
    # get iterations per depth
    for num in vbc_df['depth']:
        it_per_depth[num] += 1

    return it_per_depth

def main(args, interactive=False):
    parsed_args = parse_arguments(args)
    set_params(parsed_args)
    input = params['input'][0]
    if not os.path.exists(params['outdir']):
        os.makedirs(params['outdir'])
    if params['save']:
        vbcr.main(input,params['outdir'])
        exit()
    if os.path.isfile(input):
        if input.endswith(".vbc"):
            print("Generating single plot of instance {}".format('.'.join(input.split('.')[:-2])))
            vbc_df, treeinfo = vbcr.read(input)
            data = calc_itperdepth(vbc_df)
            plot(data,'.'.join(input.split('.')[:-2]),input.split('.')[-2], params, interactive=interactive)
        elif input.endswith("vbc.pkl"):
            vbc_df = pd.read_pickle(input)
            data = calc_itperdepth(vbc_df)
            plot(data,'.'.join(input.split('.')[:-2]),input.split('.')[-2], params, interactive=interactive)
        else:
            print("Warning: The given file is no '.vbc' or '.vbc.pkl' file.\nTerminating.")
    elif os.path.isdir(input):
        vb1 = False
        vb2 = False
        print("Generating plots from directory {}".format(input))
        for file in os.listdir(input):
            if file.endswith("vbc.pkl"):
                vbc_df = pd.read_pickle(os.path.join(input, file))
                data = calc_itperdepth(vbc_df)
                plot(data,'.'.join(file.split('.')[:-2]),file.split('.')[-2], params, interactive=interactive)
                vb1 = True
            elif file.endswith(".vbc"):
                try:
                    vbc_df, treeinfo = vbcr.read(os.path.join(input,file))
                except TypeError:
                    print("No .vbc file found. Skipping.")
                    continue
                except:
                    print("The .vbc file could not be read (e.g., due to terminated instance/file corrupted). Skipping.")
                data = calc_itperdepth(vbc_df)
                plot(data,'.'.join(file.split('.')[:-2]),file.split('.')[-2], params,interactive=interactive)
                vb2 = True
        if not (vb1 or vb2):
            print("Warning: The given directory does not contain any '.vbc' or '.vbc.pkl' files.\nTerminating.")
    else:
        print("Warning: File or directory could not be found.\nTerminating.")
        exit()


if __name__ == '__main__':
    main(sys.argv[1:])
