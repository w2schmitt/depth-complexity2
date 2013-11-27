#!/bin/bash

F=$1
TEST_DIR=$2
VALUE=$3
MAX=0
MAX_FILE=""


for file in $(ls $TEST_DIR/*/$VALUE/$F) 
do 

	LOCAL_MAX=$(awk 'END {print $1}' $file)
	
	if [ $MAX -lt $LOCAL_MAX ] 
	then
		echo $file $LOCAL_MAX
		MAX=$LOCAL_MAX
		MAX_FILE=$file
	fi
done

echo "$MAX_FILE $MAX"