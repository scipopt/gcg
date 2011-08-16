#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

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

TESTSETS=""
for i in `ls -1 $FILES | sed 's!\(.*\)check\.\([^ .]*\)\.\([^ ]*\)\.res!\2!g' | sort -u`
do
    TESTSETS="$TESTSETS $i"
done

export LC_NUMERIC=C

for i in $TESTSETS
do
    echo
    echo ====vvvv==== $i ====vvvv====
    awk -f cmpres.awk $AWKARGS texcmpfile="cmpres.$i.tex" `ls -1f $FILES | grep "$i\..*\.res"`
    echo ====^^^^==== $i ====^^^^====
done
