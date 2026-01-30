#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       *
#*                         Zuse Institute Berlin (ZIB)                       *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*
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
