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
export LANG=C

AWKARGS=""
DIR=`dirname $1`
EVALFILE=`basename $1 .eval`

OUTFILE=$DIR/$EVALFILE.out
ERRFILE=$DIR/$EVALFILE.err
RESFILE=$DIR/$EVALFILE.res
TEXFILE=$DIR/$EVALFILE.tex
PAVFILE=$DIR/$EVALFILE.pav

echo > $OUTFILE
echo > $ERRFILE
echo create overall output and error file
for i in `cat $DIR/$EVALFILE.eval` DONE
  do
  if test "$i" = "DONE"
      then
      break
  fi

  FILE=$i.out
  if test -e $FILE
      then
      cat $FILE >> $OUTFILE
  fi

  FILE=$i.err
  if test -e $FILE
      then
      cat $FILE >> $ERRFILE
  fi
done


TSTNAME=`echo $EVALFILE | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g'`

if test -f testset/$TSTNAME.solu
    then
    SOLUFILE=testset/$TSTNAME.solu
else if test -f testset/all.solu
    then
    SOLUFILE=testset/all.solu
else
    SOLUFILE=""
fi
fi
awk -f check_cbc.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" $AWKARGS $SOLUFILE $OUTFILE | tee $RESFILE


