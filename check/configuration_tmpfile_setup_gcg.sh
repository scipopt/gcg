#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2019 Operations Research, RWTH Aachen University       *
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

### resets and fills a batch file TMPFILE to run SCIP with
### sets correct limits, reads in settings, and controls
### display of the solving process

# environment variables passed as arguments
INSTANCE=$1      #  instance name to solve
GCGPATH=$2       # - path to working directory for test (usually, the check subdirectory)
TMPFILE=$3       # - the batch file to control SCIP
SETNAME=$4       # - specified basename of settings file, or 'default'
MSETNAME=$5      # - specified basename of master settings file, or 'default'
SETFILE=$6       # - instance/settings specific set file
THREADS=$7       # - the number of LP solver threads to use
FEASTOL=$8       # - feasibility tolerance, or 'default'
TIMELIMIT=$9     # - time limit for the solver
MEMLIMIT=${10}   # - memory limit for the solver
NODELIMIT=${11}  # - node limit for the solver
LPS=${12}        # - LP solver to use
DISPFREQ=${13}   # - display frequency for chronological output table
REOPT=${14}      # - true if we use reoptimization, i.e., using a difflist file instead if an instance file
STATISTICS=${15}
MODE=${16}

#args=("$@")
#for ((i=0; i < $#; i++)) {
#   echo "argument $((i+1)): ${args[$i]}"
#}

# new environment variables after running this script
# -None

#set solfile
#SOLFILE=$CLIENTTMPDIR/${USER}-tmpdir/$SOLBASENAME.sol

# reset TMPFILE
echo > $TMPFILE

if test $SETNAME != "default"
then
	echo set load $SETTINGS            >>  $TMPFILE
fi
if test $MSETNAME != "default"
then
    echo set loadmaster $MSETTINGS     >>  $TMPFILE
fi

# set non-default feasibility tolerance
if test $FEASTOL != "default"
then
    echo set numerics feastol $FEASTOL >> $TMPFILE
fi

# if permutation counter is positive add permutation seed (0 = default)
if test $p -gt 0
then
    echo set misc permutationseed $p   >> $TMPFILE
fi

# avoid solving LPs in case of LPS=none
if test "$LPS" = "none"
then
    echo set lp solvefreq -1           >> $TMPFILE
fi
echo set limits time $TIMELIMIT        >> $TMPFILE
echo set limits nodes $NODELIMIT       >> $TMPFILE
echo set limits memory $MEMLIMIT       >> $TMPFILE
echo set lp advanced threads $THREADS  >> $TMPFILE
echo set timing clocktype 1            >> $TMPFILE
echo set display verblevel 4           >> $TMPFILE
echo set display freq $DISPFREQ        >> $TMPFILE
# avoid switching to dfs - better abort with memory error
echo set memory savefac 1.0            >> $TMPFILE
echo set save $SETFILE                 >> $TMPFILE
if test "$STATISTICS" = "true"
then
    echo set visual vbcfilename $GCGPATH/results/vbc/$SHORTPROBNAME.$SETNAME.vbc >> $TMPFILE
fi

if test "$REOPT" = false
then
    # read and solve the instance
	echo read $PROB                    >> $TMPFILE

    if test $MODE = "detect"
    then
        echo presolve                      >> $TMPFILE
        echo detect                        >> $TMPFILE
        echo display statistics            >> $TMPFILE
        if test $STATISTICS = "true"
        then
            echo display additionalstatistics  >> $TMPFILE
        fi
	elif test $MODE = "miplibfeaturesoriginal"
    then
		if test ! -e results/features_original
		then
		    mkdir results/features_original
		fi
		echo set detection consclassifier consnamelevenshtein enabled FALSE    >> $TMPFILE
		echo set detection consclassifier consnamenonumbers enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevalues enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevaluesigns enabled FALSE    >> $TMPFILE
		echo set detection varclassifier scipvartype enabled FALSE    >> $TMPFILE
		echo set write miplib2017features TRUE  >> $TMPFILE
		echo set write miplib2017featurefilepath results/features_original/featurefile >> $TMPFILE
                echo detect                        >> $TMPFILE
		echo quit                          >> $TMPFILE
            elif test $MODE = "miplibfeaturespresolved"
            then
		if test ! -e results/features_presolved
		then
		    mkdir results/features_presolved
		fi
		echo set detection consclassifier consnamelevenshtein enabled FALSE    >> $TMPFILE
		echo set detection consclassifier consnamenonumbers enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevalues enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevaluesigns enabled FALSE    >> $TMPFILE
		echo set detection varclassifier scipvartype enabled FALSE    >> $TMPFILE

		echo set write miplib2017features TRUE  >> $TMPFILE
		echo set write miplib2017featurefilepath results/features_presolved/featurefile >> $TMPFILE
                echo presolve                      >> $TMPFILE
                echo detect                        >> $TMPFILE
		echo quit                          >> $TMPFILE
            elif test $MODE = "miplibfeaturesplotsoriginal"
            then
		if test ! -e $DIR/features_original
		then
		    mkdir $DIR/features_original
		fi
		if test ! -e $DIR/features_original/decs
		then
		    mkdir $DIR/features_original/decs
		fi
		if test ! -e $DIR/features_original/matrix
		then
		    mkdir $DIR/features_original/matrix
		fi
		echo set detection consclassifier consnamelevenshtein enabled FALSE    >> $TMPFILE
		echo set detection consclassifier consnamenonumbers enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevalues enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevaluesigns enabled FALSE    >> $TMPFILE
		echo set detection varclassifier scipvartype enabled FALSE    >> $TMPFILE
		echo set presolving maxrounds 0    >> $TMPFILE
		echo set write miplib2017features TRUE  >> $TMPFILE
		echo set write miplib2017plotsanddecs TRUE  >> $TMPFILE
		echo set write miplib2017shortbasefeatures TRUE  >> $TMPFILE
#		echo set display verblevel 5 >> $TMPFILE

		echo set write miplib2017featurefilepath $DIR/features_original/featurefile >> $TMPFILE
		echo set write miplib2017matrixfilepath $DIR/features_original/matrix >> $TMPFILE
		echo set write miplib2017decompfilepath $DIR/features_original/decs >> $TMPFILE
		echo set visual colorscheme 1 >> $TMPFILE
                echo detect                        >> $TMPFILE
		echo write problem  $DIR/features_original/decs/$NAME.dec     >> $TMPFILE
		echo write problem  $DIR/features_original/decs/$NAME.gp      >> $TMPFILE
		echo quit                          >> $TMPFILE
            elif test $MODE = "miplibfeaturesplotspresolved"
            then
		if test ! -e $DIR/features_presolved
		then
		    mkdir $DIR/features_presolved
		fi
		if test ! -e $DIR/features_presolved/decs
		then
		    mkdir $DIR/features_presolved/decs
		fi
		if test ! -e $DIR/features_presolved/matrix
		then
		    mkdir $DIR/features_presolved/matrix
		fi

		echo set detection consclassifier consnamelevenshtein enabled FALSE    >> $TMPFILE
		echo set detection consclassifier consnamenonumbers enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevalues enabled FALSE    >> $TMPFILE
		echo set detection varclassifier objectivevaluesigns enabled FALSE    >> $TMPFILE
		echo set detection varclassifier scipvartype enabled FALSE    >> $TMPFILE

		echo set write miplib2017features TRUE  >> $TMPFILE
		echo set write miplib2017plotsanddecs TRUE  >> $TMPFILE
		echo set write miplib2017shortbasefeatures TRUE  >> $TMPFILE

		echo set write miplib2017featurefilepath $DIR/features_presolved/featurefile >> $TMPFILE
		echo set write miplib2017matrixfilepath $DIR/features_presolved/matrix >> $TMPFILE
		echo set write miplib2017decompfilepath $DIR/features_presolved/decs >> $TMPFILE
		echo set visual colorscheme 1 >> $TMPFILE
		echo presolve                      >> $TMPFILE
                echo detect                        >> $TMPFILE
		echo write trans  $DIR/features_presolved/decs/$NAME.dec     >> $TMPFILE
		echo write trans  $DIR/features_presolved/decs/$NAME.gp      >> $TMPFILE
		echo quit                          >> $TMPFILE

	    elif test $MODE = "detectionstatistics"
	    then
		echo set detection allowclassifier enabled TRUE >> $TMPFILE
		echo presolve                      >> $TMPFILE
		echo detect                        >> $TMPFILE
                echo display detectionst           >> $TMPFILE
	    elif test $MODE = "checkexistence"
	    then
		echo set presolving maxrounds 0    >> $TMPFILE
		echo presolve                      >> $TMPFILE
		echo detect                        >> $TMPFILE
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
                        echo read $BLKFILE         >> $TMPFILE
                    fi
            elif test $MODE = "bip"
            then
                echo presolve                      >> $TMPFILE
                echo write prob bip\/$NAME-dec.bip >> $TMPFILE
                echo display statistics            >> $TMPFILE
                if test $STATISTICS = "true"
                then
                    echo display additionalstatistics  >> $TMPFILE
                fi


	elif test $MODE = "detectall"
	then
		echo presolve                      >> $TMPFILE
		echo detect                        >> $TMPFILE
		mkdir -p $GCGPATH/decs/$TSTNAME.$SETNAME
		mkdir -p $GCGPATH/images/$TSTNAME.$SETNAME
		echo write all $GCGPATH/decs\/$TSTNAME.$SETNAME dec >> $TMPFILE
		echo write all $GCGPATH/images\/$TSTNAME.$SETNAME gp >> $TMPFILE
	else
		if test -f $DECFILE -a $MODE = "readdec"
		then
		    if test -f $DECFILE
		    then
			    BLKFILE=$DECFILE
		    fi
		    if test -f $BLKFILE
		    then
			    presol=`grep -A1 PRESOLVE $BLKFILE`
		        # if we find a presolving file
			    if test $? = 0
			    then
                    # look if its in there
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
      if test "$STATISTICS" = "true"
      then
         echo display additionalstatistics  >> $TMPFILE
      fi
#       echo display solution              >> $TMPFILE
		echo checksol                      >> $TMPFILE
	fi
else
    # read the difflist file
	cat $GCGPATH/$INSTANCE                >> $TMPFILE
fi

echo quit                              >> $TMPFILE
