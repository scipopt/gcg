#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       *
#*                         Zuse Institute Berlin (ZIB)                       *
#*                                                                           *
#* This program is free software; you can redistribute it and/or             *
#* modify it under the terms of the GNU Lesser General Public License        *
#* as published by the Free Software Foundation; either version 3            *
#* of the License, or (at your option) any later version.                    *
#*                                                                           *
#* This program is distributed in the hope that it will be useful,           *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
#* GNU Lesser General Public License for more details.                       *
#*                                                                           *
#* You should have received a copy of the GNU Lesser General Public License  *
#* along with this program; if not, write to the Free Software               *
#* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#
#@file    evalheurs.sh
#@brief   Script to generate statistics for GCG Heuristics
#@author  Christian Puchert
#
# Usage: ./evalheurs.sh [ masterheurs=... ] [ origheurs=...  ] [ SETFILE(S) ] OUTFILE
#
# This script can be called on a GCG output file and will produce a .heu file
# containing statistics about the running times / calls / found sol's of the specified heuristics.
# The heuristics to be considered are specified by e. g. masterheurs=greedycolsel origheurs=gcgrounding,gcgshifting
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
	     AWKARGS="$AWKARGS -v $i"
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
