#!/usr/bin/env python3
import plotter_general as plotter
import sys
import os.path
import numpy as np
import argparse

# Checks if the bounds are either a float or the string "max"
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    if value == "max":
         return True
    else:
        return False

filters=[]
alltimes=["TOTAL TIME","PRICING TIME","PRICING SOLVER TIME","HEURISTICS TIME ORIG","HEURISTICS TIME MASTER","DETECTION TIME","PRESOLVING TIME","ORIGINAL LP TIME","RMP LP TIME","FARKAS TIME","CUTS TIME ORIG","CUTS TIME MASTER","ROOT NODE TIME","READING TIME","COPYING TIME"]

def parse_arguments(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-A', '--all',
                        action="store_true",
                        help='create all visualizations')

    parser.add_argument('-B', '--bar',
                        action="store_true",
                        help='create barchart')

    parser.add_argument('-C', '--compare',
                        action="store_true",
                        help='create comparing barchart')

    parser.add_argument('-s', '--single',
                        action="store_true",
                        help='create a pie and bar chart plot of a single instance/plot the sum of all instances in the outfile')

    parser.add_argument('-G', '--grouped-bar',
                        action="store_true",
                        help='create barchart, grouped by CPU times')

    parser.add_argument('-P', '--plot',
                        action="store_true",
                        help='create simple plot')

    parser.add_argument('-t', '--times', nargs='*',
    					default=["PRICING TIME","RMP LP TIME"],
    					help='times to be plotted. "OTHERS" time will be added automatically. The first argument is the one to be sorted by. If you want to plot all possible times, just type ALL.')

    parser.add_argument('--show-all-times', action="store_true",
                        help="shows all times that are possible")

    parser.add_argument('--buckets', type=int,
    					default=-1,
    					help='amount of buckets (resolution of grouped_bar plot)')

    parser.add_argument('--linearbucketing', action='store_true',
    					help='use linearbucketing instead of "equal amount"-bucketing')

    parser.add_argument('-load', '--loadpickle', action='store_true',
                        default=False,
                        help='Load a pickle, do not parse outfile.')

    parser.add_argument('-save', '--savepickle', action='store_true',
                        default=False,
                        help='Save a pickle, do not generate visualizations.')

    parser.add_argument('filename', nargs='+',
                        help='.out-files or .pkl-files to create plots with')

    args = parser.parse_args(args)
    if not args.savepickle and (not args.all and not args.bar and not args.grouped_bar and not args.plot and not args.compare and not args.single):
        print('All plots are disabled. Exiting script.')
        exit()

    # Check if in the times list, there are filters given
    times = []
    i = 0
    while i < len(args.times):
        if i+2 < len(args.times) and (isinstance(args.times[i], str) and isfloat(args.times[i+1]) and isfloat(args.times[i+2])):
            if args.times[i+2] == "max":
                filters.append([args.times[i],float(args.times[i+1]),np.inf])
            else:
                filters.append([args.times[i],float(args.times[i+1]),float(args.times[i+2])])
            times.append(args.times[i])
            i += 3
        else:
            times.append(args.times[i])
            i += 1
    args.times = times

    return args


def parse_outfiles_to_pickle(outfiles, save=False, path=""):
    print("Parsing .out file(s)...")
    import parser_general as pklparser
    return [pklparser.main(outfiles,save=save, path=path)]

def open_pickle(pkls):
    # Read pickle
    for file in pkls:
        if os.path.isfile(file):
            try:
                plotter.read_pickle(file)
            except:
                print("File could not be read.")
                exit()
        else:
            print("File not found.")
            exit()

def set_times(times):
    if "ALL" in times:
        times = alltimes
    # Add "other" time to list
    elif not "TOTAL TIME" in times:
        times.insert(0,"TOTAL TIME")
    return times

def set_filters(times):
    # Filter out all NaNs
    for time in times:
        plotter.add_filter(time, 0, np.inf)
    # Apply user's filters
    for filter in filters:
        plotter.add_filter(filter[0],filter[1],filter[2])
    # Some predefined filters
    #plotter.add_filter("STATUS", 1, 2)
    #plotter.add_filter("TOTAL TIME", 2, np.inf)

def plot(times, buckets, all, bar, grouped_bar, plot, compare,single, files, out, linearbucketing=False):
    if not compare:
        try:
            file=files[0].split("/")[-1].split(".")[1] + "." + files[0].split("/")[-1].split(".")[-3]
        except:
            file = "plot"
    else:
        file = files.copy()
        for i in range(len(files)):
            file[i] = files[i].split("/")[-1].split(".")[1] + "." + files[i].split("/")[-1].split(".")[-3]
    if bar or all:
    	plotter.time_distribution(times, buckets, title="TimeDistribution", type='bar', outdir=out, filename=file)
    if grouped_bar or all:
    	plotter.time_distribution(times, buckets, title="TimeDistribution", type='grouped_bar', outdir=out, filename=file, linearbucketing=linearbucketing)
    if plot or all:
    	plotter.time_distribution(times, buckets, title="TimeDistribution", type='plot', outdir=out, filename=file)
    if compare and len(files) == 1:
        print("Information: --compare given with one file. Using --single mode.")
        compare = False
        single = True
    if compare:
        print("Information: Generating comparison plots.")
        npkl = plotter.sum_pickles(files)
        plotter.time_distribution(times, npkl, title="TimeDistribution", type='compare', outdir=out, filename=files)
        plotter.time_distribution(times, npkl, title="TimeDistribution", type='comparepie', outdir=out, filename=files)
    if single:
        npkl = plotter.sum_pickles(files)
        plotter.time_distribution(times, npkl, title="TimeDistribution", type='pie', outdir=out, filename=file)
        plotter.time_distribution(times, npkl, title="TimeDistribution", type='bar', outdir=out, filename=file)

def main():
    args = sys.argv[1:]
    args = parse_arguments(args)

    # If show all times, only show times
    if args.show_all_times:
        print("\nCurrently supported times are:")
        for time in alltimes:
            print(" " + time)
        exit()

    # If saving was enabled, save to pickle and terminate
    if args.savepickle:
        parse_outfiles_to_pickle(args.filename, save=True, path=args.outdir)
        print("Done with saving as {}.general.pkl.".format(os.path.join(args.outdir, args.filename[0].split('/')[-1])))
        exit()

    # If user forgot -p for pickle file, set it on true if some file is a pickle
    for file in args.filename:
        if file.endswith(".pkl") and not args.loadpickle:
            print("Interpreting your input as pickle file")
            args.loadpickle = True

    # Parse pickles, else set pickle
    if not args.loadpickle: #thus, the out has to be parsed to pickle first
        pkl = parse_outfiles_to_pickle(args.filename, save=False)[0]
        plotter.append_open_pickle(pkl)
    else:
        pkl = args.filename
        open_pickle(pkl)

    # User warnings
    if args.buckets != -1 and not args.grouped_bar and not args.all:
        print("Information: You set a bucket amount but don't plot a grouped_bar plot.")

    # Set times, filters, pathes
    args.times = set_times(args.times)
    if not args.compare: set_filters(args.times)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Start plotting
    plot(args.times, args.buckets, args.all, args.bar, args.grouped_bar, args.plot, args.compare, args.single, args.filename,args.outdir,args.linearbucketing)

if __name__ == '__main__':
    main()
