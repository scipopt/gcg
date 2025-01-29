#!/usr/bin/env bash
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
# @author Gerald Gamrath
#
# This is what you have to do to get DIP running:
# 1) download dip, e.g. "svn checkout https://projects.coin-or.org/svn/Dip/trunk dip"
# 2) cd dip
# 3) ./configure --with-blas=BUILD --with-lapack=BUILD
#       --with-cplex-incdir=<path to your CPLEX include files, e.g., /xyz/cplex124/cplex/include/ilcplex>
#       --with-cplex-lib=<path to the CPLEX library, e.g., /xyz/cplex124/cplex/lib/x86-64_sles10_4.1/static_pic>
#       --disable-cplex-libcheck --enable-gnu-packages
# 4) make
# 5) make install
# 6) cd Dip/examples/MILPBlock
# 7) make to create the decomp_milpblock executable
# 8) create a softlink named "dip" in the check folder of gcg pointing to the
#    decomp_milpblock executable
# 9) make TEST=xyz TIME=abc ... testdip in the gcg main directory
#
# DIP can only read in mps files and the block structure definition files for DIP have
# to be created. To do this, run scripts/createBlockFiles.py xyz.test with xyz being the
# testset you want to use. This converts all instances of the testset to the mps format,
# writes them to an instance in the same folder, .lp(.gz) replaced by .mps, and also
# creates the block file in that folder, ending with .mps.block.
# The parameter file for DIP is stored at check/dip.parm, note that BlockFile and Instance
# are set by a command line overwrite.
#

TSTNAME=$1
DIPBIN=$2
SETNAME=$3
BINNAME=$DIPBIN.$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
CONTINUE=${10}

echo "timelimit: ${TIMELIMIT}"

echo "binname: ${BINNAME}"

if test ! -e results
then
    mkdir results
fi
if test ! -e settings
then
    mkdir settings
fi

GCGPATH=`pwd`

OUTFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.parm

# set the path to the settings file
if test $SETNAME != "default"
then
    SETTINGS=../settings/$SETNAME.parm
else
    SETTINGS=../settings/dip.parm
fi

if test "$CONTINUE" = "true"
then
    MVORCP=cp
else
    MVORCP=mv
fi

DATEINT=`date +"%s"`
if test -e $OUTFILE
then
    $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
fi
if test -e $ERRFILE
then
    $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
fi

if test "$CONTINUE" = "true"
then
    LASTPROB=`awk -f getlastprob.awk $OUTFILE`
    echo Continuing benchmark. Last solved instance: $LASTPROB
    echo "" >> $OUTFILE
    echo "----- Continuing from here. Last solved: $LASTPROB -----" >> $OUTFILE
    echo "" >> $OUTFILE
else
    LASTPROB=""
fi

uname -a >>$OUTFILE
uname -a >>$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

# we add 10% to the hard time limit and additional 10 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 10\` + \`expr $TIMELIMIT / 10\``

# since bash counts cpu time we need the time limit for each thread
HARDTIMELIMIT=`expr $HARDTIMELIMIT \* $THREADS`

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`

echo "hard time limit: $HARDTIMELIMIT s" >>$OUTFILE
echo "hard mem limit: $HARDMEMLIMIT k" >>$OUTFILE

for i in `cat testset/$TSTNAME.test`
do
    if test "$LASTPROB" = ""
    then
        DIR=`dirname $i`
        NAME=`basename $i .gz`
        NAME=`basename $NAME .mps`
        NAME=`basename $NAME .lp`
        DECFILE=$DIR/$NAME.dec
        LASTPROB=""
        if test -f $i
        then
            rm -f $SETFILE
            echo @01 $i ===========
            echo @01 $i ===========                 >> $ERRFILE
            echo @05 SETTINGS: $SETNAME

	    #i=`echo "$i" | sed 's/.gz//'`
	    #i=`echo "$i" | sed 's/.lp/.mps/'`

	    echo $i

	    cp $SETTINGS $TMPFILE

	    # change the time limit in the param file
	    sed -i "s,\$TIMELIMIT,$TIMELIMIT," $TMPFILE

	    # change the time limit in the param file
	    sed -i "s,\$THREADS,$THREADS," $TMPFILE

	    # change the time limit in the param file
	    sed -i "s,\$NODELIMIT,$NODELIMIT," $TMPFILE

	    # change the time limit in the param file
	    sed -i "s,\$GCGPATH,$GCGPATH," $TMPFILE

	    # change the instance in the param file
	    sed -i "s,\$INSTANCE,$i," $TMPFILE

	    # change the block file in the param file
	    sed -i "s,\$BLOCK,$DECFILE," $TMPFILE

            cp $TMPFILE $SETFILE
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            date +"@03 %s"
	    echo bash -c "ulimit -t $HARDTIMELIMIT; ulimit -v $HARDMEMLIMIT; ulimit -f 1000000; $DIPBIN --param $SETFILE"
            bash -c "ulimit -t $HARDTIMELIMIT; ulimit -v $HARDMEMLIMIT; ulimit -f 1000000; $DIPBIN --param $SETFILE" 2>>$ERRFILE
            date +"@04 %s"
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            echo =ready=
        else
            echo @02 FILE NOT FOUND: $i ===========
            echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
        fi
    else
        echo skipping $i
        if test "$LASTPROB" = "$i"
        then
            LASTPROB=""
        fi
    fi
done | tee -a $OUTFILE

rm -f $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

./evalcheck_dip.sh $OUTFILE
