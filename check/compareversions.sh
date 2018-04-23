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
if [[ $GLOBALFLAGS = *"TEST="* ]]; then
	# First cut the testset name from a copy of the global flags
	TESTNAME="${GLOBALFLAGS}"
	TESTNAME=${TESTNAME#*TEST=}
	TESTNAME=${TESTNAME%% *}

	# Store testset (Here is your problem if there are copy issues, name is hardcoded.)
	mkdir -p testset
	cd testset
	cp "$TESTNAME".test "$TESTNAME"_comparecopy.test
	cd ..

	# Replace testset name in global flags by copy
	GLOBALFLAGS=${GLOBALFLAGS//"$TESTNAME"/"$TESTNAME"_comparecopy}
fi

# If a global settings file was specified, check whether it exists.
SETTINGSNAME="${GLOBALFLAGS}"
SETTINGSNAME=${SETTINGSNAME#*SETTINGS=}
SETTINGSNAME=${SETTINGSNAME%% *}

if [ ! -f ../settings/${SETTINGSNAME}.set ]; then
    echo "Warning: Global setting file ${SETTINGSNAME}.set not found! GCG will use default settings instead."
fi

# TODO Add new folder for all files generated in the current run (we are currently in check directory)
RESDIR="results/compareversions$(date '+%d-%m-%Y_%H-%M')"
mkdir -p $RESDIR

# Script is in check, so switch to gcg main folder
cd ..
index=0
while [ $index -lt $nversions ]
do
	index=$((index + 1))

	# If a settings file for this version was specified, check whether it exists.
	SETTINGSNAME="${ADDFLAGS[$index]}"
	SETTINGSNAME=${SETTINGSNAME#*SETTINGS=}
	SETTINGSNAME=${SETTINGSNAME%% *}

	if [ ! -f ../settings/${SETTINGSNAME}.set ]; then
	    echo "Warning: Additional setting file ${SETTINGSNAME}.set not found! GCG will use default settings instead."
	fi

	# get version
	git submodule foreach --recursive git clean -f
	git submodule foreach git reset --hard
	git checkout "${VERSION[$index]}"
	git submodule init
	git submodule sync
	git submodule update
	make soplex
	make scip
	make deps ${GLOBALFLAGS} ${ADDFLAGS[$index]}
	make -j ${GLOBALFLAGS} ${ADDFLAGS[$index]}

	# 'make test' assumes there are links from gcg/check/instances to the instance folders.
	# if your testset instances are not found make sure the links are set
	make test ${GLOBALFLAGS} ${ADDFLAGS[$index]}

	# change name of output files: sort by last modified and take the first one
	cd check/results
	
	# establish VERSIONNAME as a name safe to use as a filename
	VERSIONNAME="${VERSION[$index]//\/}"		# remove slashs
	VERSIONNAME="${VERSIONNAME//.}"			# remove points
	# add a number to the version if there exists a file with the same name
	numbers='^[0-9]+$'
	firsttry=true	# bool to make sure a name ending on a number is not changed
	while [ -f "../${RESDIR}/${VERSIONNAME}.res" ]
	do
		echo ${VERSIONNAME} exists
		lastchar="${VERSIONNAME: -1}"
		if [ firsttry ] && [[ $lastchar =~ $numbers ]] ; then
			VERSIONNAME=${VERSIONNAME}1
			firsttry=false
		else
			if [[ $lastchar =~ $numbers ]] ; then
				VERSIONNAME=${VERSIONNAME::-1}
				lastchar=$((lastchar + 1))
				VERSIONNAME=${VERSIONNAME}${lastchar}
			else
				VERSIONNAME=${VERSIONNAME}1
			fi
		fi
	done
	
	# move the resulting files to the result directory of this comparison run
	OLDout=$(find . -type f -name "*.out"  -printf "%p\n" | sort -n | head -n 1)
	mv "$OLDout" "../${RESDIR}/${VERSIONNAME}.out"

	OLDres=$(find . -type f -name "*.res"  -printf "%p\n" | sort -n | head -n 1)
	mv "$OLDres" "../${RESDIR}/${VERSIONNAME}.res"

	OLDerr=$(find . -type f -name "*.err"  -printf "%p\n" | sort -n | head -n 1)
	mv "$OLDerr" "../${RESDIR}/${VERSIONNAME}.err"

	OLDtex=$(find . -type f -name "*.tex"  -printf "%p\n" | sort -n | head -n 1)
	mv "$OLDtex" "../${RESDIR}/${VERSIONNAME}.tex"

	OLDpav=$(find . -type f -name "*.pav"  -printf "%p\n" | sort -n | head -n 1)
	mv "$OLDpav" "../${RESDIR}/${VERSIONNAME}.pav"

	OLDset=$(find . -type f -name "*.set"  -printf "%p\n" | sort -n | head -n 1)
	mv "$OLDset" "../${RESDIR}/${VERSIONNAME}.set"

	# go back to the main folder to check out next version correctly
	cd ../..
done

# Return to branch the script was called on
git checkout "${CURRENTBRANCH}"

# Remove copy of global testset
cd check/testset
rm "$TESTNAME"_comparecopy.test
cd ..

# 3) do sth with the output:
echo "Start plotting..."

mkdir -p pickles
chmod +x parseres.py
chmod +x plotcomparedres.py

# parse the res files to a readable format
shopt -s nullglob
for i in ${RESDIR}/*.res; do
	./parseres.py $i ${RESDIR}
done

# plot the results
./plotcomparedres.py ${RESDIR}

echo "Finished."
echo "The results and plots can be found in $RESDIR." 

# termination
exit 0
