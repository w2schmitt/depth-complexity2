#!/bin/bash

MODEL_DIR=$1
STEP_LIMIT=$2
RESOLUTION=$3
TEST_DIR="Tests"

EXEC_NAME="depthcomplexity3d_offline"

for m in $MODEL_DIR/*
do
    NAME="$(basename ${m})"
    echo "Current model is ${NAME}"
    CURR_DIR=${TEST_DIR}/${NAME}_${RESOLUTION}
    DSTEPS=$STEP_LIMIT
    ITH=1
    #BEST=0
    ./${EXEC_NAME} -f ${MODEL_DIR}/${NAME} -dsteps ${DSTEPS} -res ${RESOLUTION} -fh "hist.txt" -fr "bestRays.off" -it 1 > info.txt

    #until [ $DSTEPS -gt $STEP_LIMIT ]; do
    #    echo "Current DSTEPS is ${DSTEPS}"
    #   
    #    ITH=$((`awk 'END { print NR }' hist.txt`-2))
    #    if [ ${BEST[i]} -lt $ITH ]; then
    #            BEST[i]=$ITH
    #    fi
    #    ITH=$((ITH-(ITH/5)))
    mkdir -p  ${CURR_DIR}/${DSTEPS}
    mv hist.txt bestRays.off info.txt ${CURR_DIR}/${DSTEPS}
    #    let DSTEPS+=10
    #done
    #echo "BEST: ${BEST}" > ${CURR_DIR}/comp.txt
done
