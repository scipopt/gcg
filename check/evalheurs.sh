#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: $
#
#@file    evalheurs.sh
#@brief   Script to generate statistics for GCG Heuristics
#@author  Christian Puchert
#
# Usage: ./evalheurs.sh [ -v heuristics=... ] [ SETFILE(S) ] OUTFILE
#
# This script can be called on a GCG output file and will produce a .heu file
# containing statistics about the running times / calls / found sol's of all included heuristics.
# If the user provides the names of some heuristics as an input argument, only those will
# be considered in the statistics (e. g. -v heuristics=gcgrounding,gcgshifting will produce
# statistics only for those two heuristics).
# Additionally, if a settings file is provided, the script will recognize if a heuristic is
# switched off and ignore it.

export LANG=C

AWKARGS=""
FILES=""
SETFILES=""
for i in $@
do
    if test ! -e $i
    then
	     AWKARGS="$AWKARGS $i"
    else
    if test `echo $i | gawk -F . '{print $NF}'` != "set" # temporary workaround; try do achieve this without awk
    then
	     FILES="$FILES $i"
	 else
	     SETFILES="$SETFILES $i"
    fi
    fi
done

for i in $FILES
do
    NAME=`basename $i .out`
    DIR=`dirname $i`
    HEUFILE=$DIR/$NAME.heu
    OUTFILE=$DIR/$NAME.out

    TSTNAME=`echo $NAME | sed 's/checkheurs.\([a-zA-Z0-9_]*\).*/\1/g'`

    if test -f $TSTNAME.test
    then
	TESTFILE=$TSTNAME.test
    else
	TESTFILE=""
    fi

    if test -f $TSTNAME.solu
    then
	SOLUFILE=$TSTNAME.solu
    else if test -f all.solu
    then
	SOLUFILE=all.solu
    else
        SOLUFILE=""
    fi
    fi

    gawk -f heurs.awk $AWKARGS $TESTFILE $SOLUFILE $SETFILES $OUTFILE | tee $HEUFILE
done
