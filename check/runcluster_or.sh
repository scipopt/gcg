#!/usr/bin/env zsh
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
#
# @author Martin Bergner
# @author Benedikt Meier

CLIENTTMPDIR=$1
CONTINUE=$2
BINNAME=$3
TLIMIT=$4
EVALFILE=$5
JOBFILE=$6
HARDMEMLIMIT=$7
ULIMITMEMLIMIT=$8
SOLVERPATH=$9
Process=${10}
Process=$(( $Process + 1 ))

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR
then
    echo "Skipping test since the path for the tmp-dir does not exist."
    exit
fi

ulimit -v $ULIMITMEMLIMIT
ulimit -m $ULIMITMEMLIMIT

export ILOG_LICENSE_FILE=$HOME/access.ilm

BASENAME=`awk "NR==$Process" $EVALFILE`
BASENAME=`basename $BASENAME`
OUTFILE=$CLIENTTMPDIR/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/$BASENAME.err
TMPFILE=$SOLVERPATH/results/$BASENAME.tmp
FILENAME=`awk "NR==$Process" $JOBFILE`
#if test ! -f $CLIENTTMPDIR/hmetis
#then
#    cp $HOME/bin/hmetis $CLIENTTMPDIR/
#fi
if test -f $CLIENTTMPDIR/hmetis
then
     rm $CLIENTTMPDIR/hmetis
fi

# check if we want to continue
if test "$CONTINUE" = "true"
then
    if test -e $SOLVERPATH/results/$BASENAME.out
    then
        exit
    fi
fi

export PATH=$PATH:$HOME/bin/
cd $CLIENTTMPDIR
uname -a                            > $OUTFILE
uname -a                            > $ERRFILE
echo "@01 $FILENAME ==========="    >> $OUTFILE
echo "@01 $FILENAME ==========="    >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE
$SOLVERPATH/../$BINNAME < $TMPFILE   >> $OUTFILE 2>>$ERRFILE
date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo "=ready="                      >> $OUTFILE

mv $OUTFILE $SOLVERPATH/results/$BASENAME.out
mv $ERRFILE $SOLVERPATH/results/$BASENAME.err

rm -f $TMPFILE