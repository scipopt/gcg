#!/bin/bash

#####################################################################################
# README
#
# This script will run different versions of GCG using the test script for comparison.
# To this end, it checks out different branches from the git, compiles, tests, and reports.
#
# Comparing different GCG versions (that might have different compile conventions) is basically one huge workaround, so beware!
# Some older versions might require manual linking of libraries (especially tag < v2*). If the test runs overnight it might stop until someone manually presses Enter.
# You have been warned.
#
# Call this script with arguments: "global flags" "gitversion1" "flags for gitversion1" "gitversion2" "flags for gitversion2" "gitversion3" ...
# e.g. ./compareversions.sh "TEST=mytestset SETTINGS=mysettings -j" "master" "" "mybranch" "LPS=cpx"
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
ORIGINALPARAMS=$@
ninputs=$#
GLOBALFLAGS=$1
ninputs=$((ninputs - 1))
shift
nversions=0

# output a little help message
if [[ $ninputs < 3 ]]; then
cat <<EOF
Call this script with arguments: "global flags" "gitversion1" "flags for gitversion1" "gitversion2" "flags for gitversion2" "gitversion3" ...
e.g. ./compareversions.sh "TEST=mytestset SETTINGS=mysettings -j" "master" "" "mybranch" "LPS=cpx"

Put "" around your arguments if they include spaces.
Add "" after a git version for no additional flags

Useful git branches to compare are, among others:
EOF
git branch -a |grep origin/release-v |sort |cut -d / -f 3
exit 0
fi

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

TESTNAME=""
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
	# Store solu file if one exists
	if [ -f "$TESTNAME".solu ]; then
		cp "$TESTNAME".solu "$TESTNAME"_comparecopy.solu
	fi
	cd ..

	# Replace testset name in global flags by copy
	GLOBALFLAGS=${GLOBALFLAGS//"$TESTNAME"/"$TESTNAME"_comparecopy}
fi

# If a global settings file was specified, check whether it exists.
if [[ $GLOBALFLAGS = *"SETTINGS="* ]]; then
	SETTINGSNAME="${GLOBALFLAGS}"
	SETTINGSNAME=${SETTINGSNAME#*SETTINGS=}
	SETTINGSNAME=${SETTINGSNAME%% *}

	if [ ! -f ../settings/${SETTINGSNAME}.set ]; then
		echo "Warning: Global setting file ${SETTINGSNAME}.set not found! GCG will use default settings instead."
	fi
fi

# Add new folder for all files generated in the current run (we are currently in check directory)
RESDIR="results/compareversions$(date '+%d-%m-%Y_%H-%M')_${TESTNAME}"
mkdir -p $RESDIR

# Add readme with parameters
if [ ! -e $RESDIR/readme.txt ]; then
	echo "This directory contains the results of the GCG version comparison run with parameters:" > $RESDIR/readme.txt
	echo "$ORIGINALPARAMS" >> $RESDIR/readme.txt
	if [ ! -z $TESTNAME ]; then
		echo "Testset $TESTNAME" >> $RESDIR/readme.txt
	fi
fi

# Script is in check, so switch to gcg main folder
cd ..

# make sure all relevant GCG linkings exist
cd lib/
mkdir -p include
mkdir -p static
ln -s scip-git/ scip
ln -s bliss/ blissinc
ln -s bliss/libbliss.a libbliss.a
ln -s /opt/cplex/cplex/include/ilcplex/ cpxinc
ln -s /opt/cplex/cplex/lib/x86-64_linux/static_pic/libcplex.a libcplex.linux.x86_64.gnu.a
ln -s /opt/cplex/cplex/lib/x86-64_linux/static_pic/libcplex.a libcplex.linux.x86_64.gnu.so
ln -s googletest/include/ gtest
ln -s googletest/build/libgtest.a libgtest.a
ln -s ../bliss-git/ include/bliss
ln -s ../bliss-git/libbliss.a static/libbliss.a
ln -s ../googletest-git/include/gtest include/gtest
ln -s ../googletest-git/build/libgtest.a static/libgtest.a
cd .. # main folder

# add links to instance folders
cd check/instances
for i in /opt/instances/* ; do
  if [ -d "$i" ]; then
    ln -s $i $(basename "$i")
  fi
done
cd ../..

index=0
while [ $index -lt $nversions ]
do
	index=$((index + 1))

	# If a settings file for this version was specified, check whether it exists.
	if [[ ${ADDFLAGS[$index]} = *"SETTINGS="* ]]; then
		SETTINGSNAME="${ADDFLAGS[$index]}"
		SETTINGSNAME=${SETTINGSNAME#*SETTINGS=}
		SETTINGSNAME=${SETTINGSNAME%% *}
	fi

	if [ ! -f ../settings/${SETTINGSNAME}.set ]; then
		echo "Warning: Additional setting file ${SETTINGSNAME}.set not found! GCG will use default settings instead."
	fi

	# get submodules
	git submodule foreach --recursive git clean -f
	git submodule foreach git reset --hard
	git checkout "${VERSION[$index]}"
	git submodule init
	git submodule sync
	git submodule update

	# set scip links depending on scip tags
	cd lib/scip-git
	SCIPTAG=$(git describe --tags)
	TAGFRAG=${SCIPTAG#v} # remove v in front of tag
	first=${TAGFRAG:0:1}
	oldsoplexlinks=false
	zimpllinks=false

	# the linking differences appear for < v310-6000-[...]
	if [ ${first} -gt 2 ]; then
		if [ $first = 3 ]; then
			TAGFRAG=${TAGFRAG#${first}}
			if [ ${TAGFRAG:0:1} -lt 1 ]; then
				oldsoplexlinks=true
				zimpllinks=true
			else
				TAGFRAG=${TAGFRAG#*-}
				TAGFRAG=${TAGFRAG%-*}
				if [ "$TAGFRAG" -lt 6000 ]; then
					oldsoplexlinks=true
				fi
				if [ "$TAGFRAG" -lt 10000 ]; then
					zimpllinks=true
				fi
			fi
		fi
	else
		oldsoplexlinks=true
		zimpllinks=true
	fi

	# function to remove old linking, assumes to be in gcg/lib/scip-git
	function removeoldlinking {
		rm -f lib/spxinc
		rm -f lib/cpxinc
		rm -f lib/libsoplex.linux.x86_64.gnu.opt.a
		rm -f lib/libsoplex.linux.x86_64.gnu.opt.so
		rm -f lib/libcplex.linux.x86_64.gnu.a
		rm -f lib/zimplinc/zimpl
		rm -f lib/libzimpl.linux.x86_64.gnu.opt.a
		rm -f lib/libzimpl.linux.x86_64.gnu.opt.so
		rm -f lib/include/spxinc
		rm -f lib/include/cpxinc
		rm -f lib/static/libsoplex.linux.x86_64.gnu.opt.a
		rm -f lib/static/libcplex.linux.x86_64.gnu.a
		rm -f lib/include/zimplinc/zimpl
		rm -f lib/static/libsoplex.linux.x86_64.gnu.opt.a libsoplex.linux.x86_64.gnu.opt.a
	}

	# make links
	removeoldlinking
	mkdir -p lib/
	cd lib/
	if $oldsoplexlinks ; then
		ln -s ../../soplex-git/src/ spxinc
		ln -s ../../soplex-git/lib/libsoplex.linux.x86_64.gnu.opt.a libsoplex.linux.x86_64.gnu.opt.a
		ln -s ../../soplex-git/lib/libsoplex.linux.x86_64.gnu.opt.a libsoplex.linux.x86_64.gnu.opt.so
		ln -s /opt/cplex1271/cplex/include/ilcplex/ cpxinc
		ln -s /opt/cplex1271/cplex/lib/x86-64_linux/static_pic/libcplex.a libcplex.linux.x86_64.gnu.a
		ln -s /opt/scipoptsuite-3.0.0/zimpl-3.3.0/lib/libzimpl.linux.x86_64.gnu.opt.a libzimpl.linux.x86_64.gnu.opt.a
		ln -s /opt/scipoptsuite-3.0.0/zimpl-3.3.0/lib/libzimpl.linux.x86_64.gnu.opt.a libzimpl.linux.x86_64.gnu.opt.so
		mkdir -p zimplinc
		cd zimplinc/
		ln -s /opt/scipoptsuite-3.0.0/zimpl-3.3.0/src/ zimpl
		cd ..
	else
		mkdir -p include/
		cd include/
		ln -s ../../../soplex-git/src/ spxinc
		ln -s /opt/cplex1271/cplex/include/ilcplex/ cpxinc
		cd ..
		mkdir -p static/
		cd static/
		ln -s ../../../soplex-git/lib/libsoplex.linux.x86_64.gnu.opt.a libsoplex.linux.x86_64.gnu.opt.a
		ln -s /opt/cplex1271/cplex/lib/x86-64_linux/static_pic/libcplex.a libcplex.linux.x86_64.gnu.a
		cd ..
		if $zimpllinks ; then
			mkdir -p include/zimplinc/
			cd include/zimplinc/
			ln -s /opt/scipoptsuite-3.0.0/zimpl-3.3.0/src/ zimpl
			cd ../../static/
			ln -s /opt/scipoptsuite-3.0.0/zimpl-3.3.0/lib/libzimpl.linux.x86_64.gnu.opt.a libzimpl.linux.x86_64.gnu.opt.a
			cd ..
		fi
	fi
	cd ../../.. # exit to gcg folder (was in gcg/lib/scip-git/lib/ before)

	# Get the versions of the components
	GCGTAG=$(git describe --tags)
	cd lib/scip-git/
	SCIPTAG=$(git describe --tags)
	cd ../soplex-git/
	SOPLEXTAG=$(git describe --tags)
	cd ../../check/
	echo "GCG ${VERSION[$index]} ${index} ${GCGTAG}" >> $RESDIR/readme.txt
	echo "SCIP ${VERSION[$index]} ${index} ${SCIPTAG}" >> $RESDIR/readme.txt
	echo "SoPlex ${VERSION[$index]} ${index} ${SOPLEXTAG}" >> $RESDIR/readme.txt
	cd ..

	# building
	make soplex ${GLOBALFLAGS} ${ADDFLAGS[$index]}	# in some older versions not sufficient
	cd lib/soplex-git
	make ${GLOBALFLAGS} ${ADDFLAGS[$index]}
	cd ../..
	make scip ${GLOBALFLAGS} ${ADDFLAGS[$index]}
	make deps ${GLOBALFLAGS} ${ADDFLAGS[$index]}
	make -j ${GLOBALFLAGS} ${ADDFLAGS[$index]}

	# 'make test' assumes there are links from gcg/check/instances to the instance folders.
	# if your testset instances are not found make sure the links are set!
	make test ${GLOBALFLAGS} ${ADDFLAGS[$index]}

	# change name of output files: sort by last modified and take the first one
	mkdir -p check/results
	cd check/results
	
	# establish VERSIONNAME as a name safe to use as a filename
	VERSIONNAME="${VERSION[$index]//\/}"		# remove slashs
	VERSIONNAME="${VERSIONNAME//.}"			# remove points

	# if the version name is empty add a name
	if [ -z "$VERSIONNAME" ] ; then
		VERSIONNAME="head"
	fi

	# if the file name exist add a number until the name is unique
	lastchar=""
	if [ -f "../${RESDIR}/${VERSIONNAME}.res" ] ; then
		echo ${VERSIONNAME} exists
		VERSIONNAME=${VERSIONNAME}1
		while [ -f "../${RESDIR}/${VERSIONNAME}.res" ]
		do
			echo ${VERSIONNAME} exists
			# as we added a '1' before the loop we can safely increment the last character
			lastchar="${VERSIONNAME: -1}"
			VERSIONNAME=${VERSIONNAME::-1}
			# note: lastchar is a string, not a single char. (e.g. 9+1=10)
			lastchar=$((lastchar + 1))
			VERSIONNAME=${VERSIONNAME}${lastchar}
		done
	fi
	
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
if [ -f "$TESTNAME"_comparecopy.solu ]; then
	rm "$TESTNAME"_comparecopy.solu
fi
cd ..

# Recover SCIP linkings (current dir is gcg/check/)
cd ../lib/scip-git/
removeoldlinking
mkdir -p lib/include/
mkdir -p lib/static/
cd lib/include/
ln -sf ../../../soplex-git/src/ spxinc
cd ../static/
ln -sf ../../../soplex-git/lib/libsoplex.linux.x86_64.gnu.opt.a libsoplex.linux.x86_64.gnu.opt.a
cd ../../../../check/

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
