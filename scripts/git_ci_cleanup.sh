#!/bin/bash
# This script automatically removes all non-running CI runs.

FOLDER=$1
PRIVATE_TOKEN=$2

echo "Cleaning up folder ${FOLDER}."
# get all existing folders
IDS=$(ls ~/$FOLDER)
for id in $IDS; do
  # access API and ask for status
  STATUS=`echo $(curl -s --header "PRIVATE-TOKEN: ${PRIVATE_TOKEN}" "https://git.or.rwth-aachen.de/api/v4/projects/1/jobs/${id}") | sed -e "s/,/\n/g" | grep "status" -m1 | cut -d\" -f4`
  if [[ "$STATUS" == "failed" ]] || [[ "$STATUS" == "success" ]] || [[ "$STATUS" == "canceled" ]]; then
    # delete faield/succeeded runs
    echo " Removing CI job $id with status $STATUS."
    rm -rf $FOLDER/$id
  else
    # keep all others
    echo " Keeping CI job $id with status $STATUS."
  fi
done
