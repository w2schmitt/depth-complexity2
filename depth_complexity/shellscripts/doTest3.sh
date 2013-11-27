#!/bin/bash

if [ $# -lt 3 ]
then
    echo "<program> <MODELS_DIR> <STEPS> <RESOLUTION>"
    exit
fi


MODEL_DIR=$1
STEP_LIMIT=$2
RESOLUTION=$3
TEST_DIR="Tests"

EXEC_NAME="depthcomplexity3d_offline"

for m in $MODEL_DIR/*
do
    NAME="$(basename ${m})"
    NAME_NO_EXT=${NAME%%.*}
    echo "Current model is ${NAME}"
    CURR_DIR=${TEST_DIR}/${NAME}_${RESOLUTION}
    DSTEPS=$STEP_LIMIT
    ITH=1
    #BEST=0
    ./${EXEC_NAME} -f ${MODEL_DIR}/${NAME} -dsteps ${DSTEPS} -res ${RESOLUTION} -fb ${NAME_NO_EXT} -f2dh "hist2D" -fth "thick_hist.js" -fh "hist.js" -fr "bestRays.off" -it 1 > info.txt
    #./${EXEC_NAME} -f ${MODEL_DIR}/${NAME} -dsteps ${DSTEPS} -res ${RESOLUTION} -fb ${NAME_NO_EXT} -f2dh "hist2D" -fr "bestRays.off" -it 1 > info.txt

    #until [ $DSTEPS -gt $STEP_LIMIT ]; do
    #    echo "Current DSTEPS is ${DSTEPS}"
    #   
    #    ITH=$((`awk 'END { print NR }' hist.txt`-2))
    #    if [ ${BEST[i]} -lt $ITH ]; then
    #            BEST[i]=$ITH
    #    fi
    #    ITH=$((ITH-(ITH/5)))
    mkdir -p  ${CURR_DIR}/${DSTEPS}
    mv hist.js thick_hist.js hist2D.js bestRays.off info.txt ${CURR_DIR}/${DSTEPS}
    mv $m analysed/$NAME
    #    let DSTEPS+=10
    #done
    #echo "BEST: ${BEST}" > ${CURR_DIR}/comp.txt
done
