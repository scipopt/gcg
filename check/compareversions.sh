#!/bin/bash

# 1) Get parameters
echo ""
echo "This script will run different versions of GCG using the test script for comparison."
echo ""
echo "Please call this script with all parameters that both GCG versions should be run with:"
echo "e.g. 'compareversions TEST=mytestset SETTINGS=mysettings LPS=cpx'"
echo "Additional parameters for a single branch are adjustable in the following."
echo ""

if [ -z $1 ]
then
	echo "No parameters given, using defaults so far."
else
	ninputs=$#
	PARAMS=""
	while ((ninputs > 0))
	do
		PARAMS=("${PARAMS[*]}" "$1")
		ninputs=$((ninputs - 1))
		shift
	done
	echo "The given parameters are:"
	echo "${PARAMS[*]}"
fi

echo ""
echo "Enter the first git branch/hash/tag:"
read VERSION1
echo "Enter additional parameters for this branch (Press Enter if none):"
read ADDPARAMS1
if [ -z "$ADDPARAMS1" ]
then
	ADDPARAMS1=""
fi

echo ""
echo "Enter the second git branch/hash/tag:"
read VERSION2
echo "Enter additional parameters for this branch (Press Enter if none):"
read ADDPARAMS2
if [ -z "$ADDPARAMS2" ]
then
	ADDPARAMS2=""
fi

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

git checkout "${CURRENTBRANCH}"


# 3) do sth with the output: TODO

# termination
exit 0
