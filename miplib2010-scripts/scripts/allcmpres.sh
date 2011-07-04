#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: allcmpres.sh,v 1.3 2011/05/30 18:47:21 bzfheinz Exp $

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
for i in `ls -1 $FILES | sed 's!\(.*\)results/\([^ .]*\)\.\([^ ]*\)\.res!\2!g' | sort -u`
do
    TESTSETS="$TESTSETS $i"
done

echo $TESTSETS

export LC_NUMERIC=C

for i in $TESTSETS
do
    RESFILE="results/$i.res"

    echo Testset: $i > $RESFILE
    echo
    echo ====vvvv==== $i ====vvvv====
    awk -f scripts/cmpres.awk $AWKARGS `ls -1 $FILES | grep "$i\..*\.res"` | tee  -a $RESFILE
    echo ====^^^^==== $i ====^^^^====
done
