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

TMP="temporary_genPlot_toRem.gp"
DATA_TMP="temp.data"


for t in 0 1
do
    COLUMN=$(($t+1));
    DC_MAX=$(awk 'FNR==v1 {print $2}' v1=$COLUMN Tests/${DIR}/comp.txt)
    #echo $DC_MAX
#You can set the size of the plot in this line (inches):
    echo "set terminal pdf size 10, 6" > ${TMP}
    echo "set output \"${DIR}_${TERM[t]}.pdf\"" >> ${TMP}
    echo "set style fill solid 1" >> ${TMP}
    echo "set boxwidth 0.8 relative" >> ${TMP}
    echo "set xtics 1" >> ${TMP}
    echo "set ytics 1" >> ${TMP}
    echo "set yrange [0:${DC_MAX}]" >> ${TMP}
    #echo "set ydata int" >> ${TMP}

    #if [ $# -eq 2 ]; then
        echo "set key off" >> ${TMP}
    #else
    # echo "set key outside" >> ${TMP}
        #echo "set key below samplen 1" >> ${TMP}
    #fi
    echo "set title \"${DIR}_${TERM[t]} MDC Behavior\"" >> ${TMP}
    echo "set xlabel \"Discretization Steps\" " >> ${TMP}
    echo "set ylabel \"Maximum Depth Complexity\"" >> ${TMP}
    #echo "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)" >> ${TMP}
    #echo "set xrange [*:*] reverse" >> ${TMP}
    DIV_CAR="plot"
    
    #COUNT=0
    MAX=$( ls Tests/${DIR}/${TYPE[t]} | sort -nr | head -n 1 )
    MIN=$( ls Tests/${DIR}/${TYPE[t]} | sort -n | head -n 1 )
    let MAX=MAX-MIN
    echo "" > ${DATA_TMP}
    for d in $( ls Tests/${DIR}/${TYPE[t]} | sort -n );
    do
     #You can set the color function here:
        #VALUE[1]=$(echo "obase=16; scale=0; (0)" | bc -l)
        #VALUE[2]=$(echo "obase=16; scale=0; (255)" | bc -l)
        #VALUE[3]="00"
        
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


RGB=(0 0 0)
#COR="#${VALUE[1]}${VALUE[2]}${VALUE[3]}"
        #echo "${COR}"
        #echo -n "${DIV_CAR} \"Tests/${DIR}/${TYPE[t]}/$d/hist.txt\" using 1:(log10(\$2+1)+1) every ::1 w lines lt rgb \"${COR}\" title '$d'" >> ${TMP}
DIV_CAR=","
echo -n "${d} " >> ${DATA_TMP}
MDC=`tail -1 Tests/${DIR}/${TYPE[t]}/$d/hist.txt | awk '{ print $1 }'`
echo -n "${MDC} " >> ${DATA_TMP}
INTENSITY=$(echo "scale=7; (${MDC}/${DC_MAX})" | bc -l)
COMMAND="evaluate_color.rb ${INTENSITY}"
#echo $INTENSITY
T=$(ruby ${COMMAND})
IND=0
for w in $T
do
  RGB[IND]=$w
  let IND=IND+1
done

COR=$(echo "scale=0; (${RGB[0]}*65536+${RGB[1]}*256+${RGB[2]})" | bc -l)
#VALUE[2]=$(echo "obase=16; scale=0; (${RGB[1]})" | bc -l)
#VALUE[3]=$(echo "obase=16; scale=0; (${RGB[2]})" | bc -l)

#VALUE[1]=$(echo "obase=16; scale=0; (${RGB[0]})" | bc -l)
#VALUE[2]=$(echo "obase=16; scale=0; (${RGB[1]})" | bc -l)
#VALUE[3]=$(echo "obase=16; scale=0; (${RGB[2]})" | bc -l)
#for ind in 1 2 3
#        do
#if [ ${#VALUE[ind]} == 1 ];
#         then
#VALUE[ind]="0${VALUE[ind]}"
#         fi
#done

#echo ${RGB[0]},${RGB[1]},${RGB[2]}
#COR="#${VALUE[1]}${VALUE[2]}${VALUE[3]}"
#echo $COR
#echo $COR
echo "${COR} " >> ${DATA_TMP}
#echo "${RGB[1]} " >> ${DATA_TMP}
#echo "${RGB[2]}" >> ${DATA_TMP}
#echo $DATA_TMP
done 
echo "plot \"${DATA_TMP}\" using 1:2:(\$3) with boxes lc rgb variable notitle" >> ${TMP}

    echo "set output \"trash.pdf\"" >> ${TMP}
    
    gnuplot ${TMP}
    rm ${TMP}
    rm ${DATA_TMP}
    rm trash.pdf
done
