#!/bin/bash
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

export LANG=C

AWKARGS=""
FILES=""
for i in $@
do
    if test ! -e $i
    then
        AWKARGS="$AWKARGS $i"
    else
        FILES="$FILES $i"
    fi
done

export LC_NUMERIC=C

for i in $FILES
do
    NAME=`basename $i .out`
    DIR=`dirname $i`
    OUTFILE=$DIR/$NAME.out
    RESFILE=$DIR/$NAME.res
    TEXFILE=$DIR/$NAME.tex
    PAVFILE=$DIR/$NAME.pav

    TSTNAME=`echo $NAME | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g'`
    echo "$TSTNAME"
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

    awk -f check.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" $AWKARGS $TESTFILE $SOLUFILE $OUTFILE | tee $RESFILE
done
