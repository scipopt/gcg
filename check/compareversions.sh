#!/bin/bash

#####################################################################################
# README
#
# This script will run different versions of GCG using the test script for comparison.
#
# Call this script with arguments: "global flags" "gitversion1" "flags for gitversion1" "gitversion2" "flags for gitversion2" "gitversion3" ...
# e.g. ./compareversions "TEST=mytestset SETTINGS=mysettings -j" "master" "" "mybranch" "LPS=cpx"
#
# Put "" around your arguments if they include spaces.
# Add "" after a git version for no additional flags
#
# In detail:
# 1st argument: all flags to be set for all git versions 
# 2nd argument: git hash/branch/tag of 1st GCG version
# 3rd argument: additional flags for the 1st GCG version
# 4th argument: git hash/branch/tag of 2nd GCG version
# 5th argument: additional flags for the 2nd GCG version
# ...
#
#####################################################################################


# 1) Get arguments

echo ""
echo "This script will run different versions of GCG using the test script for comparison."
echo ""

# Store arguments
ninputs=$#
GLOBALFLAGS=$1
ninputs=$((ninputs - 1))
shift
nversions=0
while ((ninputs > 0))
do
	nversions=$((nversions + 1))
	VERSION[$nversions]="$1"
	ninputs=$((ninputs - 1))
	shift
	if ((ninputs != 0))
	then
		ADDFLAGS[$nversions]="$1"
		ninputs=$((ninputs - 1))
		shift
	fi
done


# 2) check out the version(s), compile, run with corresponding parameter(s)

# Store current branch to return in the end
CURRENTBRANCH=$(git symbolic-ref -q HEAD)
CURRENTBRANCH=${CURRENTBRANCH##refs/heads/}
CURRENTBRANCH=${CURRENTBRANCH:-HEAD}

# Script is in check, so switch to gcg main folder
cd ..
index=0
while [ $index -lt $nversions ]
do
	index=$((index + 1))
	# get version
	git checkout "${VERSION[$index]}"
	git submodule init
	git submodule sync
	git submodule update
	make soplex
	make scip
	make deps ${GLOBALFLAGS} ${ADDFLAGS[$index]}
	make -j ${GLOBALFLAGS} ${ADDFLAGS[$index]}

	# run testset
	make test ${GLOBALFLAGS} ${ADDFLAGS[$index]}

	# change name of output files: sort by last modified and take the first one
	cd check/results
	OLDout=$(find . -type f -name "*.out"  -printf "%p\n" | sort -n | head -n 1)
	echo "$OLDout"
	mv "$OLDout" "version${index}.out"

	OLDres=$(find . -type f -name "*.res"  -printf "%p\n" | sort -n | head -n 1)
	echo "$OLDres"
	mv "$OLDres" "version${index}.res"

	OLDerr=$(find . -type f -name "*.err"  -printf "%p\n" | sort -n | head -n 1)
	echo "$OLDerr"
	mv "$OLDerr" "version${index}.err"

	OLDtex=$(find . -type f -name "*.tex"  -printf "%p\n" | sort -n | head -n 1)
	echo "$OLDtex"
	mv "$OLDtex" "version${index}.tex"

	OLDpav=$(find . -type f -name "*.pav"  -printf "%p\n" | sort -n | head -n 1)
	echo "$OLDpav"
	mv "$OLDpav" "version${index}.pav"

	OLDset=$(find . -type f -name "*.set"  -printf "%p\n" | sort -n | head -n 1)
	echo "$OLDset"
	mv "$OLDset" "version${index}.set"
	
	# parse the res file to a readable format for later use
	cd ..
	makedir -p pickles
	chmod +x parseres.py
	./parseres.py results/version${index}.res

	# go back to the main folder to check out next version correctly
	cd ..
	
done

# Return to branch the script was called on
git checkout "${CURRENTBRANCH}"


# 3) do sth with the output: TODO


# termination
exit 0
