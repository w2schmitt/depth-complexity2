#!/bin/bash

DIR=$1

#Limits for NORMAL
LOWER_LIMIT[0]=1
UPPER_LIMIT[0]=2000
#Limits for RANDOM
LOWER_LIMIT[1]=1
UPPER_LIMIT[1]=2000

TYPE[0]="Normal"
TYPE[1]="Random"

TERM[0]="u"
TERM[1]="r"

DATA_TMP[0]="$Data/{DIR}_Normal.data"
DATA_TMP[1]="$Data/{DIR}_Random.data"



DC_MAX[0]=$(awk 'FNR==1 {print $2}' Tests/${DIR}/comp.txt)
DC_MAX[1]=$(awk 'FNR==2 {print $2}' Tests/${DIR}/comp.txt)




for t in 0 1
do
    echo "DCM: ${DC_MAX[t]} " > ${DATA_TMP[t]}
    #echo "RAND_DCM: ${DC_MAX[1]} " >> ${DATA_TMP[t]}
    echo -n "STEPS " >> ${DATA_TMP[t]}
    echo "MDC " >> ${DATA_TMP[t]}
    #echo "$DC_MAX: {DC_MAX}" >  ${DATA_TMP[t]}
    MAX=$( ls Tests/${DIR}/${TYPE[t]} | sort -nr | head -n 1 )
    MIN=$( ls Tests/${DIR}/${TYPE[t]} | sort -n | head -n 1 )
    let MAX=MAX-MIN
   
    for d in $( ls Tests/${DIR}/${TYPE[t]} | sort -n );
    do
	  
	let COUNT=COUNT+1
	if [ $d -lt ${LOWER_LIMIT[t]} ]; then continue; fi
	if [ $d -gt ${UPPER_LIMIT[t]} ]; then break; fi
	for ind in 1 2 3
	      do
	    if [ ${#VALUE[ind]} == 1 ];
		then
	    VALUE[ind]="0${VALUE[ind]}"
	    fi
    done


echo -n "${d} " >> ${DATA_TMP[t]}
MDC=`tail -1 Tests/${DIR}/${TYPE[t]}/$d/hist.txt | awk '{ print $1 }'`
echo "${MDC}" >> ${DATA_TMP[t]}


done 

done
