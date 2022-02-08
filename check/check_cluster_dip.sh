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
#* Copyright (C) 2010-2022 Operations Research, RWTH Aachen University       *
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
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.out
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.res
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.err

TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$BINNAME.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
DISPFREQ=${10}
CONTINUE=${11}
QUEUETYPE=${12}
QUEUE=${13}
PPN=${14}
CLIENTTMPDIR=${15}
NOWAITCLUSTER=${16}
EXCLUSIVE=${17}

# check all variables defined
if [ -z ${EXCLUSIVE} ]
then
    echo Skipping test since not all variables are defined
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAME       = $SETNAME"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "THREADS       = $THREADS"
    echo "FEASTOL       = $FEASTOL"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    echo "QUEUETYPE     = $QUEUETYPE"
    echo "QUEUE         = $QUEUE"
    echo "PPN           = $PPN"
    echo "CLIENTTMPDIR  = $CLIENTTMPDIR"
    echo "NOWAITCLUSTER = $NOWAITCLUSTER"
    echo "EXCLUSIVE     = $EXCLUSIVE"
    exit 1;
fi


# get current GCG path
GCGPATH=`pwd`

if test ! -e $GCGPATH/results
then
    mkdir $GCGPATH/results
fi

# set the path to the settings file
if test $SETNAME != "default"
then
    SETTINGS=$GCGPATH/../settings/$SETNAME.parm
else
    SETTINGS=$GCGPATH/../settings/dip.parm
fi

if test ! -e $SETTINGS
then
    echo Skipping test since the settings file $SETTINGS does not exist.
    exit
fi

# check if binary exists
if test ! -e $GCGPATH/$BINNAME
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# check if the slurm blades should be used exclusively
if test "$EXCLUSIVE" = "true"
then
    EXCLUSIVE=" --exclusive"
else
    EXCLUSIVE=""
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# we add 10% to the hard memory limit and additional 100MB to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``

# check whether there is enough memory on the host system, otherwise we need to submit from the target system
if test "$QUEUETYPE" = "srun"
then
    HOSTMEM=`ulimit -m`
    if test "$HOSTMEM" != "unlimited"
    then
        if [ `expr $HARDMEMLIMIT \* 1024` -gt $HOSTMEM ]
        then
            echo "Not enough memory on host system - please submit from target system (e.g. ssh opt201)."
            exit
        fi
    fi
fi

# in case of qsub queue the memory is measured in kB and in case of srun the time needs to be formatted
if test  "$QUEUETYPE" = "qsub"
then
    HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
elif test  "$QUEUETYPE" = "srun"
then
    MYMINUTES=0
    MYHOURS=0
    MYDAYS=0

    #calculate seconds, minutes, hours and days
    MYSECONDS=`expr $HARDTIMELIMIT % 60`
    TMP=`expr $HARDTIMELIMIT / 60`
    if test "$TMP" != "0"
    then
	MYMINUTES=`expr $TMP % 60`
	TMP=`expr $TMP / 60`
	if test "$TMP" != "0"
	then
	    MYHOURS=`expr $TMP % 24`
	    MYDAYS=`expr $TMP / 24`
	fi
    fi
    #format seconds to have two characters
    if test ${MYSECONDS} -lt 10
    then
	MYSECONDS=0${MYSECONDS}
    fi
    #format minutes to have two characters
    if test ${MYMINUTES} -lt 10
    then
	MYMINUTES=0${MYMINUTES}
    fi
    #format hours to have two characters
    if test ${MYHOURS} -lt 10
    then
	MYHOURS=0${MYHOURS}
    fi
    #format HARDTIMELIMT
    if test ${MYDAYS} = "0"
    then
	HARDTIMELIMIT=${MYHOURS}:${MYMINUTES}:${MYSECONDS}
    else
	HARDTIMELIMIT=${MYDAYS}-${MYHOURS}:${MYMINUTES}:${MYSECONDS}
    fi
fi

#define clusterqueue, which might not be the QUEUE, cause this might be an alias for a bunch of QUEUEs
CLUSTERQUEUE=$QUEUE

NICE=""
ACCOUNT="mip"

if test $CLUSTERQUEUE = "dbg"
then
    CLUSTERQUEUE="mip-dbg,telecom-dbg"
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "telecom-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "mip-dbg"
then
    ACCOUNT="mip-dbg"
elif test $CLUSTERQUEUE = "opt-low"
then
    CLUSTERQUEUE="opt"
    NICE="--nice=10000"
fi

# counter to define file names for a test set uniquely
COUNT=1

EVALFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.eval
echo > $EVALFILE

# loop over testset
for i in `cat testset/$TSTNAME.test` DONE
do
    if test "$i" = "DONE"
    then
	break
    fi

    # increase the index for the instance tried to solve, even if the filename does not exist
    COUNT=`expr $COUNT + 1`

    # check if problem instance exists
    if test -f $GCGPATH/$i
    then

        # the cluster queue has an upper bound of 2000 jobs; if this limit is
        # reached the submitted jobs are dumped; to avoid that we check the total
        # load of the cluster and wait until it is save (total load not more than
        # 1600 jobs) to submit the next job.
	if test "$NOWAITCLUSTER" != "1"
	then
	    if test  "$QUEUETYPE" != "qsub"
	    then
		echo "waitcluster does not work on slurm cluster"
	    fi
	    ./waitcluster.sh 1600 $QUEUE 200
	fi

	SHORTPROBNAME=`basename $i .gz`
	SHORTPROBNAME=`basename $SHORTPROBNAME .mps`
	SHORTPROBNAME=`basename $SHORTPROBNAME .lp`
	SHORTPROBNAME=`basename $SHORTPROBNAME .opb`

	FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME
	BASENAME=$GCGPATH/results/$FILENAME

	TMPFILE=$BASENAME.tmp
	SETFILE=$BASENAME.set

	echo $BASENAME >> $EVALFILE

            # in case we want to continue we check if the job was already performed
	if test "$CONTINUE" != "false"
	then
	    if test -e results/$FILENAME.out
	    then
		echo skipping file $i due to existing output file $FILENAME.out
		continue
	    fi
	fi

	i=`echo "$i" | sed 's/.gz//'`
	i=`echo "$i" | sed 's/.lp/.mps/'`

	echo $i

	cp $SETTINGS $TMPFILE

	# change the time limit in the param file
	sed -i "s,\$TIMELIMIT,$TIMELIMIT," $TMPFILE

	# change the instance in the param file
	sed -i "s,\$INSTANCE,$i," $TMPFILE

	# change the block file in the param file
	sed -i "s,\$BLOCK,$i.block," $TMPFILE

        cp $TMPFILE $SETFILE

        # additional environment variables needed by runcluster.sh
	export SOLVERPATH=$GCGPATH
	export BINNAME=$BINNAME
	export BASENAME=$FILENAME
	export FILENAME=$i
	export CLIENTTMPDIR=$CLIENTTMPDIR

            # check queue type
	if test  "$QUEUETYPE" = "srun"
	then
	    sbatch --job-name=GCG$SHORTPROBNAME --mem=$HARDMEMLIMIT -p $CLUSTERQUEUE -A $ACCOUNT $NICE --time=${HARDTIMELIMIT} ${EXCLUSIVE} --output=/dev/null runcluster_dip.sh
	elif test  "$QUEUETYPE" = "bsub"
	then
	    cp runcluster_aachen_dip.sh runcluster_tmp.sh
	    TLIMIT=`expr $HARDTIMELIMIT / 60`
	    ULIMITMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
	    sed -i 's,\$CLIENTTMPDIR,$TMP,' runcluster_tmp.sh
	    sed -i "s,\$BASENAME,$BASENAME," runcluster_tmp.sh
	    sed -i "s,\$BINNAME,$BINNAME," runcluster_tmp.sh
	    sed -i "s,\$FILENAME,$FILENAME," runcluster_tmp.sh
	    sed -i "s,\$TLIMIT,$TLIMIT," runcluster_tmp.sh
	    sed -i "s,\$SHORTPROBNAME,$SHORTPROBNAME," runcluster_tmp.sh
	    sed -i "s,\$HARDMEMLIMIT,$HARDMEMLIMIT," runcluster_tmp.sh
	    sed -i "s,\$ULIMITMEMLIMIT,$ULIMITMEMLIMIT," runcluster_tmp.sh
	    sed -i "s,\$SOLVERPATH,$SOLVERPATH," runcluster_tmp.sh
#	    sed -i "s,,," runcluster_tmp.sh

#	    less runcluster_aachen.sh
#	    bsub -J GCG$SHORTPROBNAME -M $HARDMEMLIMIT -q $QUEUE -W $TLIMIT -o /dev/null < runcluster_tmp.sh &
	    bsub -q $QUEUE -o error/out_$SHORTPROBNAME_%I_%J.txt < runcluster_tmp.sh &
#	    bsub -q $QUEUE -o /dev/null < runcluster_tmp.sh &
	else
            # -V to copy all environment variables
	    qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N GCG$SHORTPROBNAME -V -q $QUEUE -o /dev/null -e /dev/null runcluster_dip.sh
	fi
    else
	echo "input file "$GCGPATH/$i" not found!"
    fi
done
