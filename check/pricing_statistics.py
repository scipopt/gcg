#!/usr/bin/env python2

import sys
import os
import argparse

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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
    Parse the files and generate the plots
    :param files: List of files to be parsed
    :return:
    """
    for file in files:
        with open(file) as _file:
            node = 0
            pricing_round = 0
            stab_round = 0
            for line in _file:
                if not line.startswith("[src/pricer_gcg.cpp:"):
                    continue
                message = line.split("] statistic: ")[-1]
                if message.startswith("New pricing round at node"):
                    # continue working here...
                    continue
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
