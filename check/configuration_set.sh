#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       *
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

# configures environment variables for test runs both on the cluster and locally
# to be invoked inside a check(_cluster.sh) script
# This script cancels the process if required variables are not correctly set

# new environment variables defined by this script:
#    GCGPATH - absolute path to invocation working directory
#    SETTINGS - settings file
#    MSETTINGS - master settings file
#    SOLUFILE - .solu file for this test set, for parsing optimal solution values
#    VALGRINDCMD - the valgrind command to use

# input environment - these environment variables should be set before invoking this script
BINNAME=$1       # name of the binary
TSTNAME=$2       # name of the test set
SETNAME=$3       # setting file basename
MSETNAME=$4      # master setting file basename
TIMELIMIT=$5     # the time limit in seconds
TIMEFORMAT=$6    # the format for the time (sec or format)
MEMLIMIT=$7      # the memory limit in MB
MEMFORMAT=$8     # the format for hard memory limit (kB or MB)
VALGRIND=$9      # should valgrind be used?
STATISTICS=${10} # should statistics be printed?

# get current GCG path
GCGPATH=`pwd`

# check if binary exists
if test ! -e "$GCGPATH/../$BINNAME"
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# create results directory if it doesn't already exist
if test ! -e "$GCGPATH/results"
then
    mkdir $GCGPATH/results
fi

# create vbc directory if it doesn't already exist
if test ! -e "$GCGPATH/results/vbc" && "$STATISTICS" = "true"
then
    mkdir $GCGPATH/results/vbc
fi

# create settings directory if non-existent
if test ! -d "$GCGPATH/../settings/"
then
    echo Create directory settings
    mkdir $GCGPATH/../settings
fi

SETTINGS=$GCGPATH/../settings/$SETNAME.set

# check if the settings file exists
if test $SETNAME != "default"
then
    if test ! -e "$SETTINGS"
    then
        echo [WARNING] Settings file $SETTINGS does not exist.
#      exit
    fi
fi

MSETTINGS=$GCGPATH/../settings/$MSETNAME.set

# check if the master settings file exists
if test $MSETNAME != "default"
then
    if test ! -e "$MSETTINGS"
    then
        echo Skipping test since the settings file $MSETTINGS does not exist.
        exit
    fi
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

if test $TIMEFORMAT = "format"
then
    #format is (d-)HH:MM:SS
    TMP=`expr $HARDTIMELIMIT`
    HARDTIMELIMIT=""
    DIVISORS=(60 60 24)
    for((i=0; i<=2; i++))
    do
        printf -v HARDTIMELIMIT "%02d${HARDTIMELIMIT}" `expr ${TMP} % ${DIVISORS[i]}`
        # separate the numbers by colons except for the last (HH hours)
        if test $i -lt 2
        then
            HARDTIMELIMIT=":${HARDTIMELIMIT}"
        fi
        TMP=`expr ${TMP} / ${DIVISORS[i]}`
    done
    if test ${TMP} -gt 0
    then
        HARDTIMELIMIT=${TMP}-${HARDTIMELIMIT}
    fi
fi

# we add 10% to the hard memory limit and additional 100MB to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``

# qsub requires the memory limit to be displayed in kB
if test "$MEMFORMAT" = "kB"
then
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`
elif test "$MEMFORMAT" = "B"
then
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
fi
# check if the test run should be processed in the valgrind environment
if test "$VALGRIND" = "true"
then
    VALGRINDCMD="valgrind --log-fd=1 --leak-check=full "
else
    VALGRINDCMD=""
fi

#check if additional instance paths are given
POSSIBLEPATHS=$GCGPATH
if test -e paths.txt
then
    POSSIBLEPATHS="${POSSIBLEPATHS} `cat paths.txt`"
fi
POSSIBLEPATHS="${POSSIBLEPATHS} DONE"
echo $POSSIBLEPATHS
