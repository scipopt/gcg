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
# @author Martin Bergner
# @author Christian Puchert
# @author Gerald Gamrath

export LANG=C

REMOVE=0
AWKARGS=""
FILES=""

for i in $@
do
  if test ! -e $i
  then
      if test "$i" = "-r"
      then
          REMOVE=1
      else
          AWKARGS="$AWKARGS $i"
      fi
  else
      FILES="$FILES $i"
  fi
done

for FILE in $FILES
do

  DIR=`dirname $FILE`
  EVALFILE=`basename $FILE .eval`
  EVALFILE=`basename $EVALFILE .out`

  OUTFILE=$DIR/$EVALFILE.out
  ERRFILE=$DIR/$EVALFILE.err
  SETFILE=$DIR/$EVALFILE.set
  RESFILE=$DIR/$EVALFILE.res
  TEXFILE=$DIR/$EVALFILE.tex
  PAVFILE=$DIR/$EVALFILE.pav

  # check if the eval file exists; if this is the case construct the overall solution files
  if test -e $DIR/$EVALFILE.eval
  then
      echo > $OUTFILE
      echo > $ERRFILE
      echo create overall output and error file for $EVALFILE

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
       if test "$REMOVE" = "1"
       then
           rm -f $FILE
       fi
   else
       echo Missing $i
   fi

   FILE=$i.err
   if test -e $FILE
   then
       cat $FILE >> $ERRFILE
       if test "$REMOVE" = "1"
       then
           rm -f $FILE
       fi
   fi

   FILE=$i.set
   if test -e $FILE
   then
       cp $FILE $SETFILE
       if test "$REMOVE" = "1"
       then
           rm -f $FILE
       fi
   fi

   FILE=$i.tmp
   if test -e $FILE
        then
       if test "$REMOVE" = "1"
       then
           rm -f $FILE
       fi
   fi
      done

      if test "$REMOVE" = "1"
      then
          rm -f $DIR/$EVALFILE.eval
      fi
  fi

  # check if the out file exists
  if test -e $DIR/$EVALFILE.out
  then
      echo create results for $EVALFILE

      # detect test set
      TSTNAME=`echo $EVALFILE | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g'`

      # detect test used solver
      SOLVER=`echo $EVALFILE | sed 's/check.\([a-zA-Z0-9_-]*\).\([a-zA-Z0-9_]*\).*/\2/g'`

      echo "Testset " $TSTNAME
      echo "Solver  " $SOLVER

      if test -f testset/$TSTNAME.test
      then
          TESTFILE=testset/$TSTNAME.test
      else
          TESTFILE=""
      fi

      # look for .solu files under the name of the test, the name of the test with everything after the first "_" pt "-" stripped, and "_all"
      SOLUFILE=""
      for F in $TSTNAME ${TSTNAME%%_*} ${TSTNAME%%-*} _all
      do
          if test -f testset/${F}.solu
          then
              SOLUFILE=testset/${F}.solu
              break
          fi
      done

      if test  "$SOLVER" = "cplex"
      then
     awk -f check_cplex.awk -v "TEXFILE=$TEXFILE" $AWKARGS $SOLUFILE $OUTFILE | tee $RESFILE
      elif test  "$SOLVER" = "dip"
      then
          awk -f check_dip.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" $AWKARGS $TESTFILE $SOLUFILE $OUTFILE | tee $RESFILE
      else
          awk -f check.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" $AWKARGS $TESTFILE $SOLUFILE $OUTFILE | tee $RESFILE
      fi
  fi
done
