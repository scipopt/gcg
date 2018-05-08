#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import commands

def main(argv):
    presolved = (argv[5] == 'presolved') or (argv[5] == 'True')
    presolvedstring = ""
    olderinstancepathhelp = argv[9].split("/")
    olderinstancepath = ""
    for i in range(3, len(olderinstancepathhelp) ):
        if i == 3:
            olderinstancepath = olderinstancepath + olderinstancepathhelp[i]
        else:
            olderinstancepath = olderinstancepath + "/" + olderinstancepathhelp[i]


    if presolved:
        presolvedstring = " set presolving maxrounds 0 presolve"

    gcgcommand = argv[3] + " -c \"read " + argv[1] + " set load " + argv[4]
    gcgcommand = gcgcommand + " set write miplib2017features TRUE " + "set write miplib2017plotsanddecs TRUE " + " set write miplib2017shortbasefeatures TRUE "
    gcgcommand = gcgcommand + " set write miplib2017featurefilepath " + argv[6]
    gcgcommand = gcgcommand + " set write miplib2017matrixfilepath " + argv[7]
    gcgcommand = gcgcommand + " set write miplib2017decompfilepath " + argv[8]
    gcgcommand = gcgcommand + " set visual colorscheme 1"
    gcgcommand = gcgcommand + " change instancename " + olderinstancepath
    gcgcommand = gcgcommand +  presolvedstring
    gcgcommand = gcgcommand + " detect"

    if presolved:
        gcgcommand = gcgcommand + " write trans " + argv[8] + ".dec"
    else:
        gcgcommand = gcgcommand + " write prob " + argv[8] + ".dec"

    if presolved:
        gcgcommand = gcgcommand + " write trans " + argv[8] + ".gp"
    else:
        gcgcommand = gcgcommand + " write prob " + argv[8] + ".gp"

    gcgcommand = gcgcommand + " quit\" "

    print gcgcommand

    os.system(gcgcommand)

    pass

if __name__ == "__main__":
    if len(sys.argv) == 10:
        main(sys.argv)
    else:
        print "usage: featuresPlots.py <path_instance> <instancename> <path_gcg> <settings> <presolved> <path_features> <path_matrix_plots> <path_decomp_files> <path_orig_instance>"
