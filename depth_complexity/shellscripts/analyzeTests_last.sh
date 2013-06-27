#!/bin/bash

DIR=$1

#Limits for NORMAL
LOWER_LIMIT[0]=1
UPPER_LIMIT[0]=2000
#Limits for RANDOM
#LOWER_LIMIT[1]=1
#UPPER_LIMIT[1]=2000

#TYPE[0]="Normal"
#TYPE[1]="Random"

#TERM[0]="u"
#TERM[1]="r"

TMP="temporary_genPlot_toRem.gp"

#for t in 0 1
#do
    #You can set the size of the plot in this line (inches):
echo "set terminal pdf size 10, 6" > ${TMP}
echo "set output \"${DIR}.pdf\"" >> ${TMP}
if [ $# -eq 2 ]; then
    echo "set key off" >> ${TMP}
else
    echo "set key outside" >> ${TMP}
    #echo "set key below samplen 1" >> ${TMP}
fi
echo "set title \"${DIR} Histogram\"" >> ${TMP}
echo "set xlabel \"Depth Complexity\" " >> ${TMP}
echo "set ylabel \"Number of Rays\"" >> ${TMP}
DIV_CAR="plot"
COUNT=0
MAX=$( ls Tests/${DIR} | sort -nr | head -n 1 )
MIN=$( ls Tests/${DIR} | sort -n | head -n 1 )
let MAX=MAX-MIN
for d in $( ls Tests/${DIR} | sort -n );
do

    #You can set the color function here:
    VALUE[1]=$(echo "obase=16; scale=0; (${COUNT}*255/${MAX})" | bc -l)
    VALUE[2]=$(echo "obase=16; scale=0; (${COUNT}*255/${MAX})" | bc -l)
    VALUE[3]="00"

    let COUNT=COUNT+10
    if [ $d -lt ${LOWER_LIMIT} ]; then continue; fi
    if [ $d -gt ${UPPER_LIMIT} ]; then break; fi
    for ind in 1 2 3
    do
            if [ ${#VALUE[ind]} == 1 ];
            then
                    VALUE[ind]="0${VALUE[ind]}"
            fi
    done
    COR="#${VALUE[1]}${VALUE[2]}${VALUE[3]}"
    #echo "${COR}"
    echo $(gawk '{print $1}' Tests/${DIR}/$d/hist.txt)

    echo -n "${DIV_CAR} \"Tests/${DIR}/$d/hist.txt\" using (\$1):(\$2) every ::1 w lines lt rgb \"${COR}\" title '$d'" >> ${TMP}
        DIV_CAR=","
done
echo "" >> ${TMP}
echo "set output \"trash.pdf\"" >> ${TMP}
gnuplot ${TMP}
rm ${TMP}
rm trash.pdf
#done

