#!/bin/bash

TEST_DIR=$1
F=$2;
VALUE=$3
DIR=$4;

for file in $(ls $TEST_DIR/*/$VALUE/$F) 
do 
    NAME="$(basename ${file})"
    y=${file#*m}
    x=${y%%.*}
    NEW_NAME="m$x.js"
    echo $file
    cp $file $DIR
    mv $DIR/$NAME $DIR/$NEW_NAME
done


