#!/bin/bash

F=$1;
DIR=$2;
TEST_DIR=$3
VALUE=$4

for file in $(ls $TEST_DIR/*/$VALUE/$F) 
do 
    NAME="$(basename ${file})"
    y=${file#*m}
    x=${y%%.*}
    NEW_NAME="m$x.txt"
    cp $file $DIR
    mv $DIR/$NAME $DIR/$NEW_NAME
done


