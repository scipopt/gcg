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
#* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       *
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
MSETNAME=$4
BINID=$5
TIMELIMIT=$6
NODELIMIT=$7
MEMLIMIT=$8
THREADS=$9
FEASTOL=${10}
DISPFREQ=${11}
CONTINUE=${12}
QUEUETYPE=${13}
QUEUE=${14}
PPN=${15}
CLIENTTMPDIR=${16}
NOWAITCLUSTER=${17}
EXCLUSIVE=${18}
PERMUTE=${19}
MODE=${20}

# check all variables defined
if [ -z ${MODE} ]
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
    exit 1;
fi


# get current GCG path
GCGPATH=`pwd`

if test ! -e $GCGPATH/results
then
    mkdir $GCGPATH/results
fi

SETTINGS=$GCGPATH/../settings/$SETNAME.set

# check if the settings file exists
if test $SETNAME != "default"
then
    if test ! -e $SETTINGS
    then
        echo Skipping test since the settings file $SETTINGS does not exist.
        exit
    fi
fi

MSETTINGS=$GCGPATH/../settings/$MSETNAME.set

# check if the master settings file exists
if test $MSETNAME != "default"
then
    if test ! -e $MSETTINGS
    then
        echo Skipping test since the settings file $MSETTINGS does not exist.
        exit
    fi
fi

# check if binary exists
if test ! -e $GCGPATH/../$BINNAME
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# check if queue has been defined
if test "$QUEUE" = ""
then
    echo Skipping test since the queue name has not been defined.
    exit
fi

# check if number of nodes has been defined
if test "$PPN" = ""
then
    echo Skipping test since the number of nodes has not been defined.
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
COUNT=0

# loop over permutations
for ((p = 0; $p <= $PERMUTE; p++))
do
    # if number of permutations is positive, add postfix
    if test $PERMUTE -gt 0
    then
	EVALFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME"#p"$p.eval
	JOBFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME"#p"$p.job
    else
	EVALFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.eval
	JOBFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.job
    fi
    rm -f $EVALFILE
    touch $EVALFILE
    rm -f $JOBFILE
    touch $JOBFILE

    # loop over testset
    for i in `cat testset/$TSTNAME.test` DONE
    do
	if test "$i" = "DONE"
	then
	    break
	fi

        # increase the index for the inctance tried to solve, even if the filename does not exist
	COUNT=`expr $COUNT + 1`

	PROB=$GCGPATH/`echo $i|cut -d";" -f1`
	DECFILE=$GCGPATH/`echo $i|cut -d";" -f2`

        # check if problem instance exists
	if test -f $PROB
	then

            # the cluster queue has an upper bound of 2000 jobs; if this limit is
            # reached the submitted jobs are dumped; to avoid that we check the total
            # load of the cluster and wait until it is save (total load not more than
            # 1900 jobs) to submit the next job.
	    if test "$NOWAITCLUSTER" != "1"
	    then
		if test  "$QUEUETYPE" != "qsub"
		then
		    echo "waitcluster does not work on slurm cluster"
		fi
		./waitcluster.sh 1600 $QUEUE 200
	    fi

	    DIR=`dirname $PROB`
	    NAME=`basename $NAME .gz`
	    NAME=`basename $NAME .mps`
	    NAME=`basename $NAME .lp`

	    if test "$PROB" == "$DECFILE"
	    then
		EXT=${PROB##*.}
		if test "$EXT" = "gz"
		then
		    BLKFILE=$DIR/$NAME.blk.gz
		    DECFILE=$DIR/$NAME.dec.gz
		else
		    BLKFILE=$DIR/$NAME.blk
		    DECFILE=$DIR/$NAME.dec
		fi
	    fi


	    SHORTPROBNAME=`basename $NAME .gz`
	    SHORTPROBNAME=`basename $SHORTPROBNAME .mps`
	    SHORTPROBNAME=`basename $SHORTPROBNAME .lp`
	    SHORTPROBNAME=`basename $SHORTPROBNAME .opb`

	    # if number of permutations is positive, add postfix
	    if test $PERMUTE -gt 0
	    then
		FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME#"p"$p
	    else
		FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME
	    fi

	    BASENAME=$GCGPATH/results/$FILENAME

	    TMPFILE=$BASENAME.tmp
	    SETFILE=$BASENAME.set

	    echo $BASENAME >> $EVALFILE
	    echo $SHORTPROBNAME >> $JOBFILE

            # in case we want to continue we check if the job was already performed
	    if test "$CONTINUE" != "false"
	    then
		if test -e results/$FILENAME.out
		then
		    echo skipping file $i due to existing output file $FILENAME.out
		    continue
		fi
	    fi

	    echo > $TMPFILE
	    if test $SETNAME != "default"
	    then
		echo set load $SETTINGS            >>  $TMPFILE
	    fi
        if test $MSETNAME != "default"
        then
        echo set loadmaster $MSETTINGS     >>  $TMPFILE
        fi
	    if test $FEASTOL != "default"
	    then
		echo set numerics feastol $FEASTOL >> $TMPFILE
	    fi

	    # if permutation counter is positive add permutation seed (0 = default)
	    if test $p -gt 0
	    then
		echo set misc permutationseed $p   >> $TMPFILE
	    fi

	    echo set limits time $TIMELIMIT        >> $TMPFILE
	    echo set limits nodes $NODELIMIT       >> $TMPFILE
	    echo set limits memory $MEMLIMIT       >> $TMPFILE
	    echo set lp advanced threads $THREADS  >> $TMPFILE
	    echo set timing clocktype 1            >> $TMPFILE
	    echo set display verblevel 4           >> $TMPFILE
	    echo set display freq $DISPFREQ        >> $TMPFILE
	    echo set memory savefac 1.0            >> $TMPFILE # avoid switching to dfs - better abort with memory error
	    echo set save $SETFILE                 >> $TMPFILE
	    echo read $GCGPATH/$i                 >> $TMPFILE
#           echo presolve                         >> $TMPFILE
	    if test $MODE = "detect"
	    then
		echo presolve                      >> $TMPFILE
		echo detect                        >> $TMPFILE
		echo display statistics            >> $TMPFILE
		echo presolve                      >> $TMPFILE
	    else
		if test -f $DECFILE -a $MODE = "readdec"
		then
		    if test -f $DECFILE
		    then
			BLKFILE=$DECFILE
		    fi
		    if test -f $BLKFILE
		    then
			presol=`grep -A1 PRESOLVE $BLKFILE`
		    # if we find a presolving file
			if test $? = 0
			then
                        # look if its in there
			    if grep -xq 1 - <<EOF
$presol
EOF
			    then
				echo presolve      >> $TMPFILE
			    fi
			fi
			echo read $BLKFILE         >> $TMPFILE
		    fi
		fi
		echo optimize                      >> $TMPFILE
		echo display statistics            >> $TMPFILE
		echo display additionalstatistics  >> $TMPFILE
#               echo display solution                  >> $TMPFILE
		echo checksol                      >> $TMPFILE
	    fi
	    echo quit                              >> $TMPFILE

            # additional environment variables needed by runcluster.sh
	    export SOLVERPATH=$GCGPATH
	    export BINNAME=$BINNAME
	    export BASENAME=$FILENAME
	    export FILENAME=$i
	    export CLIENTTMPDIR=$CLIENTTMPDIR

            # check queue type
	    if test  "$QUEUETYPE" = "srun"
	    then
		sbatch --job-name=GCG$SHORTPROBNAME --mem=$HARDMEMLIMIT -p $CLUSTERQUEUE -A $ACCOUNT $NICE --time=${HARDTIMELIMIT} ${EXCLUSIVE} --output=/dev/null runcluster.sh
	    elif test  "$QUEUETYPE" = "bsub"
	    then
		cp runcluster_aachen.sh runcluster_tmp.sh
		TLIMIT=`expr $HARDTIMELIMIT / 60`
		ULIMITMEMLIMIT=`expr $HARDMEMLIMIT \* 1024000`
		sed -i 's,\$CLIENTTMPDIR,$TMP,' runcluster_tmp.sh
		sed -i "s,\$BINNAME,$BINNAME," runcluster_tmp.sh
		sed -i "s,\$TLIMIT,$TLIMIT," runcluster_tmp.sh
		sed -i "s,\$EVALFILE,$EVALFILE," runcluster_tmp.sh
		sed -i "s,\$JOBFILE,$JOBFILE," runcluster_tmp.sh
		sed -i "s,\$HARDMEMLIMIT,$HARDMEMLIMIT," runcluster_tmp.sh
		sed -i "s,\$ULIMITMEMLIMIT,$ULIMITMEMLIMIT," runcluster_tmp.sh
		sed -i "s,\$SOLVERPATH,$SOLVERPATH," runcluster_tmp.sh
#	        sed -i "s,,," runcluster_tmp.sh

#	        less runcluster_aachen.sh
#	        bsub -J GCG$SHORTPROBNAME -M $HARDMEMLIMIT -q $QUEUE -W $TLIMIT -o /dev/null < runcluster_tmp.sh &
#	        bsub -q $QUEUE -o error/out_$SHORTPROBNAME_%I_%J.txt < runcluster_tmp.sh &
#	        bsub -q $QUEUE -o /dev/null < runcluster_tmp.sh &
	    else
                # -V to copy all environment variables
		qsub -l walltime=$HARDTIMELIMIT -l mem=$HARDMEMLIMIT -l nodes=1:ppn=$PPN -N GCG$SHORTPROBNAME -V -q $QUEUE -o /dev/null -e /dev/null runcluster.sh
	    fi
	else
	    echo "input file "$GCGPATH/$i" not found!"
	fi
    done
    if test  "$QUEUETYPE" = "bsub"
    then
        bsub -J "$TSTNAME[1-$COUNT]" -q $QUEUE -o error/out_$TSTNAME_%I_%J.txt < runcluster_tmp.sh
    fi
done
