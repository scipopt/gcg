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

def parse_arguments(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir', type=str,
                        default="plots",
                        help='output directory (default: "plots")')

    parser.add_argument('-t', '--times', nargs='*',
                        default=["TOTAL TIME", "RMP LP TIME"],
                        help='times to be plotted. Use exactly two times, all others will be ignored.')

    parser.add_argument('-load', '--loadpickle', action='store_true',
                        default=False,
                        help='Load a pickle, do not parse outfile.')

    parser.add_argument('-save', '--savepickle', action='store_true',
                        default=False,
                        help='Save a pickle, do not generate visualizations.')

    parser.add_argument('filename', nargs='+',
                        help='.out-files or .pkl-files to create plots with')

    args = parser.parse_args(args)

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

    # Check if correct amount of times is given
    if len(args.times) not in (0,2):
        parser.error('Either give no values for times, or two, not {}.'.format(len(args.times)))
        exit()

    return args


def parse_outfiles_to_pickle(outfiles, save=False, path=""):
    print("Parsing .out file(s)...")
    import parser_general as pklparser
    return [pklparser.main(outfiles,save=save)]

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

def set_filters(times):
    # Filter out all NaNs
    for time in times:
        plotter.add_filter(time, 0, np.inf)
    # Some predefined filters
    plotter.add_filter("STATUS", 1, 2)
    plotter.add_filter("TOTAL TIME", 2, np.inf)
    # Set user filters
    for filter in filters:
        plotter.add_filter(filter[0],filter[1],filter[2])

def main():
    args = sys.argv[1:]
    args = parse_arguments(args)

    # If saving was enabled, save to pickle and terminate
    if args.savepickle:
        parse_outfiles_to_pickle(args.filename, save=True, path=args.outdir)
        print("Done with saving as {}.general.pkl.".format(args.filename[0]))
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

    # Set times, filters, pathes
    set_filters(args.times)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    args.filename = [args.filename[0].split('/')[-1]]
    try:
        settings = "." + args.filename[0].split(".")[15]
    except:
        settings = ".default"
    if len(args.filename) > 2:
        file = args.filename[0].split(".")[1]
    else:
        file = args.filename[0].split(".")[0]

    plotter.twinplot([args.times[0]], [args.times[1]], outdir=args.outdir, filename=file + settings)


if __name__ == '__main__':
    main()
