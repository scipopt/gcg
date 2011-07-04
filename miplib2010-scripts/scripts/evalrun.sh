#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: evalrun.sh,v 1.6 2011/04/06 20:16:06 bzfheinz Exp $

AWKARGS=""
FILES=""

# construct paths
MIPLIBPATH=`pwd`

# absolut tolerance for checking linear constraints and objective value
LINTOL=1e-4 
# absolut tolerance for checking integrality constraints 
INTTOL=1e-4 

for i in $@
do
    if test ! -e $i
    then
	AWKARGS="$AWKARGS $i"
    else
	FILES="$FILES $i"
    fi
done

for i in $FILES
do
    NAME=`basename $i .out`
    DIR=`dirname $i`
    OUTFILE=$DIR/$NAME.out
    RESFILE=$DIR/$NAME.res

    TSTNAME=`echo $NAME | sed 's/\([a-zA-Z0-9_-]*\).\([a-zA-Z0-9_-]*\).*/\1/g'`
    SOLVER=`echo $NAME | sed 's/\([a-zA-Z0-9_-]*\).\([a-zA-Z0-9_-]*\).*/\2/g'`

    SOLUFILE=$MIPLIBPATH/testsets/$TSTNAME.solu

     # check if a solution  file/link exists
    if test ! -e $SOLUFILE
    then
	echo "Warning: solution file/link <$TSTNAME.solu> does not exist in <testsets> folder; therefore, no consistency check"
	SOLUFILE=""
    fi
    
    awk -f scripts/parse.awk -f scripts/parse_$SOLVER.awk -v "LINTOL=$LINTOL" $SOLUFILE $OUTFILE | tee $RESFILE
done
