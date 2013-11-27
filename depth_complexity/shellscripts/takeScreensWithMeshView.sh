#!/bin/bash

MODELS_DIR=$1
MESH_VIEW=$2

for file in $(ls $MODELS_DIR)
do
	NAME="$(basename ${file})"
    y=${file#*m}
    x=${y%%.*}
	./${MESH_VIEW} ${x} -grab ${MODELS_DIR}/${file} 
done
