#!/bin/bash

#####################################################################################
# README
#
# This script will run different versions of GCG using the test script for comparison.
#
# Call this script with arguments: "global flags" "gitversion1" "flags for gitversion1" "gitversion2" "flags for gitversion2" "gitversion3" ...
# e.g. ./compareversions "TEST=mytestset SETTINGS=mysettings -j" "master" " " "mybranch" "LPS=cpx"
#
# Put "" around your arguments if they include spaces.
# Put " " for empty flags. (!)
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

# Sanity checks for arguments
if [ -z $1 ]
then
	echo "No arguments."
	exit 0
fi

# Store arguments
ninputs=$#
GLOBALFLAGS=$1
ninputs=$((ninputs - 1))
shift
nversions=0
while ((ninputs > 0))
do
	nversions=$((nversions + 1))
	VERSION${nversions}=$1
	ninputs=$((ninputs - 1))
	shift
	if ((ninputs != 0))
	then
		ADDFLAGS${nversions}=$1
		ninputs=$((ninputs - 1))
		shift
	fi
done

echo $nversions


# 2) TODO check out the version(s), compile, run with corresponding parameter(s)
#		if out files would get overwritten then add version/params coding to their names

# Store current branch to return in the end
CURRENTBRANCH=$(git symbolic-ref -q HEAD)
CURRENTBRANCH=${CURRENTBRANCH##refs/heads/}
CURRENTBRANCH=${CURRENTBRANCH:-HEAD}


# test version 1
cd ..					# Script is in check, so switch to gcg main folder
git checkout "${VERSION1}"
git submodule init
git submodule sync
git submodule update
make soplex
make scip
make deps ${PARAMS[*]}
make -j ${PARAMS[*]}

# run testset
make test ${PARAMS[*]}

# change name of output files
cd check/results
OLDout=$(find . -type f -name "*.out"  -printf "%p\n" | sort -n | head -n 1)
echo "$OLDout"
mv "$OLDout" "version1.out"

OLDres=$(find . -type f -name "*.res"  -printf "%p\n" | sort -n | head -n 1)
echo "$OLDres"
mv "$OLDres" "version1.res"

# Return to branch the script was called on
git checkout "${CURRENTBRANCH}"


# 3) do sth with the output: TODO

# termination
exit 0
