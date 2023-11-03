#!/usr/bin/env zsh
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       *
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
module switch intel gcc/4.6

#BSUB -J SCIP$SHORTFILENAME
#BSUB -M $HARDMEMLIMIT
#BSUB -W $TLIMIT

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR
then
    echo Skipping test since the path for the tmp-dir does not exist.
    exit
fi

ulimit -v $ULIMITMEMLIMIT
ulimit -m $ULIMITMEMLIMIT

export ILOG_LICENSE_FILE=$HOME/access.ilm

BASENAME=`awk "NR==$LSB_JOBINDEX" $EVALFILE`
BASENAME=`basename $BASENAME`
OUTFILE=$CLIENTTMPDIR/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/$BASENAME.err
TMPFILE=$SOLVERPATH/results/$BASENAME.tmp
FILENAME=`awk "NR==$LSB_JOBINDEX" $JOBFILE`
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
$SOLVERPATH/../$BINNAME < $TMPFILE  >> $OUTFILE 2>>$ERRFILE
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