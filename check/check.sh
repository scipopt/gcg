#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       *
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
# @author Martin Bergner
# @author Gerald Gamrath
# @author Christian Puchert

TSTNAME=$1
BINNAME=$2
SETNAME=$3
MSETNAME=$4
BINID=$5
TIMELIMIT=$6
NODELIMIT=$7
MEMLIMIT=$8
THREADS=$9
FEASTOL=${10}
DISPFREQ=${11}
CONTINUE=${12}
LOCK=${13}
VERSION=${14}
LPS=${15}
VALGRIND=${16}
MODE=${17}

SETDIR=../settings

if test ! -e results
then
    mkdir results
fi
if test ! -e locks
then
    mkdir locks
fi

LOCKFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.lock
RUNFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.run.$BINID
DONEFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.done

OUTFILE=results/check.$TSTNAME.$BINID.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINID.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINID.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINID.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINID.$SETNAME.set

SETTINGS=$SETDIR/$SETNAME.set
MSETTINGS=$SETDIR/$MSETNAME.set

if test "$LOCK" = "true"
then
    if test -e $DONEFILE
    then
        echo skipping test due to existing done file $DONEFILE
        exit
    fi
    if test -e $LOCKFILE
    then
        if test -e $RUNFILE
        then
            echo continuing aborted run with run file $RUNFILE
        else
            echo skipping test due to existing lock file $LOCKFILE
            exit
        fi
    fi
    date > $LOCKFILE
    date > $RUNFILE
fi

if test ! -e $OUTFILE
then
    CONTINUE=false
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

# we add 50% to the hard time limit and additional 10 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 10\` + \`expr $TIMELIMIT / 2\``

# we add 10% to the hard memory limit and additional 1000mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 1000\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`

echo "hard time limit: $HARDTIMELIMIT s" >>$OUTFILE
echo "hard mem limit: $HARDMEMLIMIT k" >>$OUTFILE

VALGRINDCMD=
if test "$VALGRIND" = "true"
then
   VALGRINDCMD="valgrind --log-fd=1 --leak-check=full"
fi

for i in `cat testset/$TSTNAME.test` DONE
do
    if test "$i" = "DONE"
    then
        date > $DONEFILE
        break
    fi

    PROB=`echo $i|cut -d";" -f1`
    DECFILE=`echo $i|cut -d";" -f2`
    DIR=`dirname $PROB`
    NAME=`basename $PROB .gz`
    NAME=`basename $NAME .mps`
    NAME=`basename $NAME .lp`

    if test "$PROB" == "$DECFILE"
    then
        EXT=${PROB##*.}
        if test "$EXT" = "gz"
        then
            BLKFILE=$DIR/$NAME.blk.gz
            DECFILE=$DIR/$NAME.dec.gz
        else
            BLKFILE=$DIR/$NAME.blk
            DECFILE=$DIR/$NAME.dec
        fi
    fi
    if test "$LASTPROB" = ""
    then
        if test -f $PROB
        then
            echo @01 $PROB ===========
            echo @01 $PROB ===========             >> $ERRFILE
            echo > $TMPFILE
            if test "$SETNAME" != "default"
            then
                echo set load $SETTINGS            >> $TMPFILE
            fi
            if test "$MSETNAME" != "default"
            then
                echo set loadmaster $MSETTINGS     >>  $TMPFILE
            fi
            if test "$FEASTOL" != "default"
            then
                echo set numerics feastol $FEASTOL >> $TMPFILE
            fi
            echo set limits time $TIMELIMIT        >> $TMPFILE
            echo set limits nodes $NODELIMIT       >> $TMPFILE
            echo set limits memory $MEMLIMIT       >> $TMPFILE
            echo set lp advanced threads $THREADS  >> $TMPFILE
            echo set timing clocktype 1            >> $TMPFILE
            echo set display verblevel 4           >> $TMPFILE
            echo set display freq $DISPFREQ        >> $TMPFILE
            echo set memory savefac 1.0            >> $TMPFILE # avoid switching to dfs - better abort with memory error
            if test "$LPS" = "none"
            then
                echo set lp solvefreq -1           >> $TMPFILE # avoid solving LPs in case of LPS=none
            fi
            echo set save $SETFILE                 >> $TMPFILE
            echo read $PROB                        >> $TMPFILE

            if test $MODE = "detect"
            then
                echo presolve                      >> $TMPFILE
                echo detect                        >> $TMPFILE
                echo display statistics            >> $TMPFILE
            elif test $MODE = "bip"
            then
                echo presolve                      >> $TMPFILE
                echo write prob bip\/$NAME-dec.bip >> $TMPFILE
                echo display statistics            >> $TMPFILE
            elif test $MODE = "detectall"
            then
                echo presolve                      >> $TMPFILE
                echo detect                        >> $TMPFILE
                mkdir -p decs/$TSTNAME.$SETNAME
                mkdir -p images/$TSTNAME.$SETNAME
                echo write all decs\/$TSTNAME.$SETNAME dec >> $TMPFILE
                echo write all images\/$TSTNAME.$SETNAME gp >> $TMPFILE
            else
                if test $MODE = "readdec"
                then
                    if test -f $DECFILE
                    then
                        BLKFILE=$DECFILE
                    fi
                    if test -f $BLKFILE
                    then
                        EXT=${BLKFILE##*.}
                        if test "$EXT" = "gz"
                        then
                            presol=`zgrep -A1 PRESOLVE $BLKFILE`
                        else
                            presol=`grep -A1 PRESOLVE $BLKFILE`
                        fi
                        echo $presol
                        # If the decomposition contains presolving information ...
                        if test $? = 0
                        then
                            # ... check if it belongs to a presolved problem
                            if grep -xq 1 - <<EOF
$presol
EOF
                            then
                                echo presolve      >> $TMPFILE
                            fi
                        fi
                        echo read $BLKFILE         >> $TMPFILE
                    fi
                fi
                echo optimize                      >> $TMPFILE
                echo display statistics            >> $TMPFILE
#               echo display additionalstatistics  >> $TMPFILE
#               echo display solution              >> $TMPFILE
                echo checksol                      >> $TMPFILE
            fi
            echo quit                              >> $TMPFILE
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            date +"@03 %s"
            bash -c " ulimit -t $HARDTIMELIMIT s; ulimit -v $HARDMEMLIMIT k; ulimit -f 200000; $VALGRINDCMD ../$BINNAME < $TMPFILE" 2>>$ERRFILE
            date +"@04 %s"
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            echo
            echo =ready=
            if test $MODE = "detectall"
            then
                mv *_*.dec decs\/
            fi
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
rm -f cipreadparsetest.cip

date >>$OUTFILE
date >>$ERRFILE

if test -e $DONEFILE
then
    ./evalcheck.sh $OUTFILE

    if test "$LOCK" = "true"
    then
        rm -f $RUNFILE
    fi
fi
