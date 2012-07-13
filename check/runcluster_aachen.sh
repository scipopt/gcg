#!/usr/bin/env zsh
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
OUTFILE=$CLIENTTMPDIR/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/$BASENAME.err
TMPFILE=$SOLVERPATH/results/$BASENAME.tmp

#if test ! -f $CLIENTTMPDIR/hmetis
#then
#    cp $HOME/bin/hmetis $CLIENTTMPDIR/
#fi
if test -f $CLIENTTMPDIR/hmetis
then
     rm $CLIENTTMPDIR/hmetis
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
#chmod g+r $ERRFILE
#chmod g+r $SCIPPATH/results/$BASENAME.out
#chmod g+r $SCIPPATH/results/$BASENAME.set
