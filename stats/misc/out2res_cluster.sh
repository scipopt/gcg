#!/usr/bin/env bash
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
# This is the evalcheck_cluster.sh script from the check folder, just without the
# arguments to also export .tex and .pav and a changed directory for the check.awk
#

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
          awk -f misc/check.awk $AWKARGS $TESTFILE $SOLUFILE $OUTFILE | tee $RESFILE
      fi
  fi
done
