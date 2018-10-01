#!/usr/bin/env bash
#set -x
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
# @author Martin Bergner
# @author Gerald Gamrath
# @author Christian Puchert
# @author Benedikt Meier
#
# Call with "make testcluster"
#
# The queue is passed via $QUEUE (possibly defined in a local makefile in scip/make/local).
#
# For each run, we can specify the number of nodes reserved for a run via $PPN. If tests runs
# with valid time measurements should be executed, this number should be chosen in such a way
# that a job is run on a single computer, i.e., in general, $PPN should equal the number of cores
# of each computer. Of course, the value depends on the specific computer/queue.
#
# To get the result files call "./evalcheck_cluster.sh
# results/check.$TSTNAME.$BINNAME.$SETNAME.eval in directory check/
# This leads to result files
#  - results/check.$TSTNAME.$BINNAME.$SETNAME.$MSETNAME.out
#  - results/check.$TSTNAME.$BINNAME.$SETNAME.$MSETNAME.res
#  - results/check.$TSTNAME.$BINNAME.$SETNAME.$MSETNAME.err

TSTNAME=$1
BINNAME=$2
SETNAME=$3
MSETNAME=$4
BINID=$5
TIMELIMIT=$6
NODELIMIT=$7
MEMLIMIT=$8
THREADS=$9
FEASTOL=${10}
LPS=${11}
DISPFREQ=${12}
CONTINUE=${13}
QUEUETYPE=${14}
QUEUE=${15}
PPN=${16}
CLIENTTMPDIR=${17}
NOWAITCLUSTER=${18}
EXCLUSIVE=${19}
PERMUTE=${20}
MODE=${21}
STATISTICS=${22}
PROJECT=${23}

# check all variables defined
if [ -z ${PROJECT} ]
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAME       = $SETNAME"
    echo "MSETNAME      = $MSETNAME"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "THREADS       = $THREADS"
    echo "FEASTOL       = $FEASTOL"
    echo "LPS           = $LPS"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    echo "QUEUETYPE     = $QUEUETYPE"
    echo "QUEUE         = $QUEUE"
    echo "PPN           = $PPN"
    echo "CLIENTTMPDIR  = $CLIENTTMPDIR"
    echo "NOWAITCLUSTER = $NOWAITCLUSTER"
    echo "EXCLUSIVE     = $EXCLUSIVE"
    echo "PERMUTE       = $PERMUTE"
    echo "MODE          = $MODE"
    echo "STATISTICS    = $STATISTICS"
    echo "PROJECT       = $PROJECT"
    exit 1;
fi

# configure cluster-related environment variables
. ./configuration_cluster.sh $QUEUE $PPN $EXCLUSIVE $QUEUETYPE

# the srun queue requires a format duration HH:MM:SS (and optionally days),
# whereas the qsub requires the memory limit in kB
if test "QUEUETYPE" != "qsub"
then
    TIMEFORMAT="format"
    MEMFORMAT="MB"
else
    TIMEFORMAT="sec"
    MEMFORMAT="B"
fi
# call routines for creating the result directory, checking for existence
# of passed settings, etc
VALGRIND="false" # change this to "true" for using valgrind
. ./configuration_set.sh $BINNAME $TSTNAME $SETNAME $MSETNAME $TIMELIMIT $TIMEFORMAT $MEMLIMIT $MEMFORMAT $VALGRIND $STATISTICS

# for RWTH computing projects
if test "$PROJECT" = "none"
then
    PROJFLAG=""
else
    PROJFLAG="-P $PROJECT"
fi

# at the first time, some files need to be initialized. set to "" after the innermost loop
# finished the first time
INIT="true"

# loop over permutations
for ((p = 0; $p <= $PERMUTE; p++))
do
    # counter to define file names for a test set uniquely
    COUNT=0

    # loop over testset
    for INSTANCE in `cat testset/$TSTNAME.test` DONE
    do
        if test "$INSTANCE" = "DONE"
        then
            break
        fi

        # increase the index for the instance tried to solve, even if the filename does not exist
	    COUNT=`expr $COUNT + 1`

	    PROB=$GCGPATH/`echo $INSTANCE|cut -d";" -f1`
	    DECFILE=$GCGPATH/`echo $INSTANCE|cut -d";" -f2`

        # check if problem instance exists
	    if [ ! -f $PROB ]
	    then
	        echo "input file "$PROB" not found!"
            continue
	    fi

        # the cluster queue has an upper bound of 2000 jobs; if this limit is
        # reached the submitted jobs are dumped; to avoid that we check the total
        # load of the cluster and wait until it is save (total load not more than
        # 1600 jobs) to submit the next job.
        if test "${NOWAITCLUSTER}" -eq "0" && ( test "$QUEUETYPE" = "qsub" || test "$QUEUETYPE" = "bsub" )
        then
            ./waitcluster.sh 1600 $QUEUE 200
        elif test "${NOWAITCLUSTER}" -eq "0"
        then
            echo "waitcluster does not work on slurm cluster"
        fi

        # infer the names of all involved files from the arguments
        . ./configuration_logfiles.sh $INIT $COUNT $INSTANCE $BINID $PERMUTE $SETNAME $MSETNAME $TSTNAME $CONTINUE $QUEUE $p
        
        # copy the basename into eval-file for invokation of evalcheck_cluster
        if test "$SKIPINSTANCE" = "true"
        then
            continue
        fi

        # call tmp file configuration for GCG
        . ./configuration_tmpfile_setup_gcg.sh $INSTANCE $GCGPATH $TMPFILE $SETNAME $MSETNAME $SETFILE $THREADS $FEASTOL $TIMELIMIT $MEMLIMIT $NODELIMIT $LPS $DISPFREQ false $CLIENTTMPDIR $STATISTICS $FILENAME

        # additional environment variables needed by runcluster.sh
	    export SOLVERPATH=$GCGPATH
	    export BINNAME=$BINNAME
	    export BASENAME=$FILENAME
	    export FILENAME=$INSTANCE
	    export CLIENTTMPDIR=$CLIENTTMPDIR

        # check queue type
        if test  "$QUEUETYPE" = "srun"
        then
            # additional environment variables needed by run.sh
            export SOLVERPATH=$GCGPATH
            export EXECNAME=$GCGPATH/../$BINNAME
            export BASENAME=$FILENAME
            export FILENAME=$INSTANCE
            export CLIENTTMPDIR
            export HARDTIMELIMIT
            export HARDMEMLIMIT
            export CHECKERPATH=$SCIPPATH/solchecker
            sbatch --job-name=GCG$SHORTPROBNAME --mem=$HARDMEMLIMIT -p $CLUSTERQUEUE -A $ACCOUNT $NICE --time=${HARDTIMELIMIT} ${EXCLUSIVE} --output=/dev/null runcluster.sh
        elif  test "$QUEUETYPE" = "qsub"
        then
            # -V to copy all environment variables
            qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N SCIP$SHORTPROBNAME -v SOLVERPATH=$SCIPPATH,EXECNAME=$SCIPPATH/../$BINNAME,BASENAME=$FILENAME,FILENAME=$i,CLIENTTMPDIR=$CLIENTTMPDIR -V -q $CLUSTERQUEUE -o /dev/null -e /dev/null runcluster.sh
        fi


        INIT="false"
    done # end for TSTNAME

    # queue type on RWTH Aachen cluster is bsub
    if test  "$QUEUETYPE" = "bsub"
    then
		cp runcluster_aachen.sh runcluster_tmp.sh
        TLIMIT=`echo $HARDTIMELIMIT | awk '{ n = split($0,a,":"); print 60*a[1]+a[2];}'`
		ULIMITMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
		sed -i 's,\$CLIENTTMPDIR,$TMP,' runcluster_tmp.sh
		sed -i "s,\$CONTINUE,$CONTINUE," runcluster_tmp.sh
		sed -i "s,\$BINNAME,$BINNAME," runcluster_tmp.sh
		sed -i "s,\$TLIMIT,$TLIMIT," runcluster_tmp.sh
		sed -i "s,\$EVALFILE,$EVALFILE," runcluster_tmp.sh
		sed -i "s,\$JOBFILE,$JOBFILE," runcluster_tmp.sh
		sed -i "s,\$HARDMEMLIMIT,$HARDMEMLIMIT," runcluster_tmp.sh
		sed -i "s,\$ULIMITMEMLIMIT,$ULIMITMEMLIMIT," runcluster_tmp.sh
		sed -i "s,\$SOLVERPATH,$SOLVERPATH," runcluster_tmp.sh

	    bsub -J "$TSTNAME[1-$COUNT]" -q $QUEUE -o error/out_$TSTNAME_%I_%J.txt $PROJFLAG < runcluster_tmp.sh
    fi

    # queue type on RWTH Aachen cluster OR is condor
    if test  "$QUEUETYPE" = "condor"
    then
        TLIMIT=`echo $HARDTIMELIMIT | awk '{ n = split($0,a,":"); print 60*a[1]+a[2];}'`
        ULIMITMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
        AGUMENTS=`echo "$CLIENTTMPDIR $CONTINUE $BINNAME $TLIMIT $EVALFILE $JOBFILE $HARDMEMLIMIT $ULIMITMEMLIMIT $SOLVERPATH"`
        condor_submit <<EOF
################################
#
################################
Universe = vanilla
# Universe = standard

should_transfer_files = IF_NEEDED
when_to_transfer_output = ON_EXIT
# when_to_transfer_output = ON_EXIT_OR_EVICT

Executable = $GCGPATH/../check/runcluster_or.sh
Arguments  = " ${AGUMENTS} \$(process)"

Log = $SOLVERPATH/results/condor_Log.\$(Cluster).\$(process).txt
Output = $SOLVERPATH/results/condor_Out.\$(Cluster).\$(process).txt
Error = $SOLVERPATH/results/condorErr.\$(Cluster).\$(process).txt

Queue $COUNT
EOF

    fi



done # end for PERMUTE
