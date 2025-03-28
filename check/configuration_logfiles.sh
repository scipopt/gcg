#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
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

### configures the right test output files such as the .eval, the .tmp and the .set
### files to run a test on.
### the invoking script should pass "init" as argument to signalize that
### files need to be reset

### environment variables declared in this script
### OUTFILE - the name of the (sticked together) output file
### ERRFILE - the name of the (sticked together) error file
### SHORTPROBNAME - the basename of $INSTANCE without file extension
### FILENAME - the basename of the local files (.out, .tmp, and .err)
### EVALFILE - evaluation file to glue single output and error files together
### SKIPINSTANCE - should the instance be skipped because it was already evaluated in a previous setting?
### BASENAME - $GCGPATH/results/$FILENAME cf. FILENAME argument
### TMPFILE  - the batch file name to pass for solver instructions
### SETFILE  - the name of the settings file to save solver settings to

### environment variables passed as arguements to this script
INIT=$1      # should log files be initialized (this overwrite or copy/move some existing log files)
COUNT=$2     # the instance count as part of the filename
INSTANCE=$3  # the name of the instance
BINID=$4     # the ID of the binary to use
PERMUTE=$5   # the number of permutations to use - 0 for no permutation
SETNAME=$6   # the name of the setting
MSETNAME=$7  # the name of the master setting
TSTNAME=$8   # the name of the testset
CONTINUE=$9  # should test continue an existing run
# optional variables
QUEUE=${10}   # the queue name
p=${11}       # the index of the current permutation - only needed if permutations are used

if test "$QUEUE" = ""
then
    QUEUE=`hostname`
fi

OUTFILE=results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME.err

# if number of permutations is positive, add postfix
if test $PERMUTE -gt 0
then
    EVALFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME"#p"$p.eval
    JOBFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME"#p"$p.job
else
    EVALFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME.eval
    JOBFILE=$GCGPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME.job
fi

if test "$INSTANCE" = "DONE"
then
    return
fi

# reset files if flag is set to 'init'
if test $INIT = "true"
then
    # reset the eval and job files
    rm -f $EVALFILE
    touch $EVALFILE
    rm -f $JOBFILE
    touch $JOBFILE

    # mv existing out and error files
    if test "$CONTINUE" = "true"
    then
        MVORCP=cp
    else
        MVORCP=mv
    fi
    DATEINT=`date +"%s"`
    for FILE in OUTFILE ERRFILE
    do
        if test -e $FILE
        then
            $MVORCP $FILE $FILE.old-$DATEINT
        fi
    done
fi

# filter all parseable file format extensions
SHORTPROBNAME=`basename $PROB .gz`
for EXTENSION in .mps .lp .opb .gms .pip .zpl .cip .fzn .osil .wbo .cnf
do
    SHORTPROBNAME=`basename $SHORTPROBNAME $EXTENSION`
done

# if no decomposition was specified, look for a .dec or .blk file in the instance directory
DIR=`dirname $PROB`
if test "$PROB" == "$DECFILE"
then
	EXT=${PROB##*.}
	if test "$EXT" = "gz"
	then
		BLKFILE=$DIR/$SHORTPROBNAME.blk.gz
		DECFILE=$DIR/$SHORTPROBNAME.dec.gz
	else
		BLKFILE=$DIR/$SHORTPROBNAME.blk
		DECFILE=$DIR/$SHORTPROBNAME.dec
	fi
fi

# if number of permutations is positive, add postfix
if test $PERMUTE -gt 0
then
    FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME#"p"$p
else
    FILENAME=$USER.$TSTNAME.$COUNT"_"$SHORTPROBNAME.$BINID.$QUEUE.$SETNAME.$MSETNAME
fi

SKIPINSTANCE="false"
# in case we want to continue we check if the job was already performed
if test "$CONTINUE" = "true" && test -e results/$FILENAME.out
then
    echo skipping file $INSTANCE due to existing output file results/$FILENAME.out
    SKIPINSTANCE="true"
fi

# configure global names TMPFILE (batch file) and SETFILE to save settings to
BASENAME=$GCGPATH/results/$FILENAME
TMPFILE=$BASENAME.tmp
SETFILE=$BASENAME.set

echo $BASENAME >> $EVALFILE
echo $SHORTPROBNAME >> $JOBFILE
