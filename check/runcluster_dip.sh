#!/usr/bin/env bash
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
# @author Gerald Gamrath

# check if tmp-path exists
if test ! -d $CLIENTTMPDIR
then
    echo Skipping test since the path for the tmp-dir does not exist.
    exit
fi

OUTFILE=$CLIENTTMPDIR/$BASENAME.out
ERRFILE=$CLIENTTMPDIR/$BASENAME.err
TMPFILE=$SOLVERPATH/results/$BASENAME.tmp

uname -a                            > $OUTFILE
uname -a                            > $ERRFILE
echo @01 $FILENAME ===========      >> $OUTFILE
echo @01 $FILENAME ===========      >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE
date +"@03 %s"                      >> $OUTFILE
$SOLVERPATH/$BINNAME --param $TMPFILE  >> $OUTFILE 2>>$ERRFILE
date +"@04 %s"                      >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
echo -----------------------------  >> $OUTFILE
date                                >> $ERRFILE
echo                                >> $OUTFILE
echo =ready=                        >> $OUTFILE

mv $OUTFILE $SOLVERPATH/results/$BASENAME.out
mv $ERRFILE $SOLVERPATH/results/$BASENAME.err

rm -f $TMPFILE
#chmod g+r $ERRFILE
#chmod g+r $SCIPPATH/results/$BASENAME.out
#chmod g+r $SCIPPATH/results/$BASENAME.set
