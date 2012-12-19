#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "Usage: $0 input.ref"
	exit 1
fi

NUMBLKS=$(cat $1 | head -n 1 | cut -d ' ' -f 1)
REFHEADER=$(cat $1 | head -n 1 | cut -d ' ' -f 2-)
for i in $(seq 1 $NUMBLKS); do
	LINE=$(cat $1 | head -n $[$i + 1] | tail -n 1)
	NUMCONS=$(echo $REFHEADER | cut -d ' ' -f $i)
	if [ -z $NUMCONS ]; then
		continue
	fi
	echo "$[i - 1] $NUMCONS"
	echo $LINE | sed 's/,/ /g'
done
