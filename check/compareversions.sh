#!/bin/bash

#####
# README: 
# The script is called with parameters of the form 
# "testset testparam1;testparam2;... githash/gitbranch1 param1;param2;... githash/gitbranch2 param1;param2;..."
# Parameters:
#	- testset:		Testset to run the different versions on
#	- testparams:		Parameters for the test script
#	- githash/gitbranchN:	GCG versions to be run (a githash or a branch)
#	- params:		parameters with which to compile & run GCG version #N
#####


# For readability reasons: 
# Even though not necessary please "declare" global variables here if they are guaranteed to exist.
TESTSET=""		# testset on which to run all versions
USERINPUT=""		# temp var for user inputs
VERSIONCOUNTER=0	# stores amount of versions to compare


# 1) Get parameters
echo ""
echo "This script will run different versions of GCG using the test script for comparison."

if [ -z $1 ]
then
	echo "No input to compare. Please call this script as"
	echo "./compareversions testset testparam1,testparam2,... githash/gitbranch1 param1,param2,... githash/gitbranch2 param1,param2,..."
	echo "Place \",,\" for empty parameters."
	exit 0
fi

# read & store params
ninputs=$#
TESTSET=$1
ninputs=$((ninputs - 1))
shift
TESTPARAMS=$1
ninputs=$((ninputs - 1))
shift
VERSION=""
PARAMS=""
nversions=0
echo $ninputs
while ((ninputs > 0))
do
	nversions=$((nversions + 1))
	VERSION=("${VERSION[*]}" "$1")
	ninputs=$((ninputs - 1))
	shift
	PARAMS=("${PARAMS[*]}" "$1")
	ninputs=$((ninputs - 1))
	shift
done
echo VERSIONCOUNTER= "$nversions", VERSION= "${VERSION[*]}", PARAMS= "${PARAMS[*]}". 

# 2) TODO check out the version(s), compile, run with corresponding parameter(s)
#		if out files would get overwritten then add version/params coding to their names

# Store current branch to return in the end
CURRENTBRANCH=$(git symbolic-ref -q HEAD)
CURRENTBRANCH=${CURRENTBRANCH##refs/heads/}
CURRENTBRANCH=${CURRENTBRANCH:-HEAD}

i=1
while (( i <= ${VERSIONCOUNTER} ))
do
	# checkout current version
	cd ..					# Script is in check, so switch to gcg main folder
	git checkout "${VERSION[${i}]}"
	git submodule init
	git submodule sync
	git submodule update
	make soplex
	make scip
	make "${COMPARAMS[${i}]}" deps
	make "${COMPARAMS[${i}]}"

	# run testset
	make test TEST="$TESTSET" SETTINGS

	# TODO change name of output files

	i=$((i + 1))
done
git checkout "${CURRENTBRANCH}"


# 3) do sth with the output: 

# termination
exit 0
