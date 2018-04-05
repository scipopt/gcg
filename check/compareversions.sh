#!/bin/bash

#####################################################################################
# README
#
# This script will run different versions of GCG using the test script for comparison.
# Comparing different GCG versions (that might have different compile conventions) is basically one huge workaround, so beware!
# Some older versions might require manual linking of libraries (especially tag < v2*). If the test runs overnight it might stop until someone manually presses Enter.
# You have been warned.
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

# If there is one test file for all versions, store the test file in case it differs in the versions
# (This is assuming that setting the test in the additional flags is intentional and test sets are supposed to differ in that case.)

# First cut the testset name from a copy of the global flags
TESTNAME="${GLOBALFLAGS}"
TESTNAME=${TESTNAME#*TEST=}
TESTNAME=${TESTNAME%% *}

# Store testset (Here is your problem if there are copy issues, name is hardcoded.)
mkdir -p testset
cd testset
cp "$TESTNAME".test "$TESTNAME"_comparecopy.test
cd ..

# TODO Replace testset name in global flags by copy
GLOBALFLAGS=${GLOBALFLAGS//"$TESTNAME"/"$TESTNAME"_comparecopy}

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

	# make sure the testsets can be found: make symbolic link (TODO for striplibn for now, there should be a more elegant version -> discuss design!)
	cd check/instances
	ln -s /opt/instances/striplibn striblibn
	cd ../..
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

	# go back to the main folder to check out next version correctly
	cd ../..
	
done

# Return to branch the script was called on
git checkout "${CURRENTBRANCH}"

# TODO Remove copy of global testset
cd check/testset
rm "$TESTNAME"_comparecopy.test
cd ..

# 3) do sth with the output: TODO
	
# parse the res files to a readable format for later use
mkdir -p pickles
chmod +x parseres.py

index=0
while [ $index -lt $nversions ]
do
	index=$((index + 1))
	./parseres.py results/version${index}.res
done

./plotcomparedres.py

# termination
exit 0
