#!/usr/bin/env bash
#
# This bash script updates all copyrights in the GCG files and posted
# those files which contain a copyright which do not have the right format
#
# You just have to run this script. There is nothing to adjust. 
# The correct year is detected through the 'date' function 
#
# Note that not all files (usually scripts) contain a copyright. A copyright is only 
# needed for those files which are part of a GCG distribution (see makedist.sh)

NEWYEAR=`date +"%Y"`
LASTYEAR=`expr $NEWYEAR - 1`

DIRECTORIES=(check doc lint src src/* stats stats/* tests)
EXTENSIONS=(sh awk h c hpp cpp html lnt py dag)
EXTRAFILES=(Makefile INSTALL check/eval make/make.project make/make.detecthost)

echo ""
echo "This script reports *all* files which have not a correct COPYRIGHT."
echo "Only files which are included in the distribution need a COPYRIGHT (see makedist.sh)"
echo ""

# collect all files
FILES=""
for DIRECTORY in ${DIRECTORIES[@]}
do
    # exclude cppad subdirectory
    if [[ "$DIRECTORY" =~ "src/cppad" ]]
    then
	continue
    fi
    for EXTENSION in ${EXTENSIONS[@]}
    do
	for FILE in $DIRECTORY/*.$EXTENSION
	do
	    if test -f $FILE
	    then
		FILES="$FILES $FILE"
	    fi
	done
    done
    for EXTRAFILE in ${EXTRAFILES[@]}
    do
	if test -f $DIRECTORY/$EXTRAFILE
	then
	    FILES="$FILES $DIRECTORY/$EXTRAFILE"
	fi
    done
done

for EXTRAFILE in ${EXTRAFILES[@]}
do 
    if test -f $EXTRAFILE
    then
	FILES="$FILES $EXTRAFILE"
    fi
done

for FILE in ${FILES[@]}
do
    if test -f $FILE
    then
	echo $FILE
	
	# check if the file has a correct old copyright 
	COUNT1=`grep -c 2010-$LASTYEAR $FILE`

	# check if the file has a correct new copyright 
	COUNT2=`grep -c 2010-$NEWYEAR $FILE`

	if [ $COUNT1 -eq 0 ] && [ $COUNT2 -gt 0 ]
	then
	    continue
	elif [ $COUNT1 -gt 0 ] && [ $COUNT2 -eq 0 ]
	then
	    echo "fixed"
	    sed -i 's!2010-'$LASTYEAR'!2010-'$NEWYEAR'!g' $FILE
	else
	    echo "COPYRIGHT ERROR ===> $FILE"
	fi
    fi
done
