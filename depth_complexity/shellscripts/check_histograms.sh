#!/bin/bash

if [ $# -lt 2 ]
then
	echo "<program> <HIST_DIR> <NUMBER_OF_MODELS> <OPT_COPY_FROM> <OPT_COPY_TO>"
	exit
fi

DIR=$1
NUM=$2
COPY_FROM=$3
COPY_TO=$4

for i in $(seq 0 $NUM) 
do 
	FILE="${DIR}/m${i}.txt"
    if [ ! -f $FILE ] 
    then
		NAME="m${i}.off"
		echo $NAME
		if [ $# -gt 3 ] 
		then
			cp $COPY_FROM/$NAME $COPY_TO/
		fi
	fi
done


