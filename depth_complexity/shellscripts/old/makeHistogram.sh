#!/bin/bash

DIR=$1

TMP="temporary_genPlot_toRem.gp"
COUNT=0
SUM=0

for FILEPATH in $DIR/*.txt
do
		FILE="${FILEPATH##*/}"
    echo "Processing ${FILE} file..."
		
    #You can set the size of the plot in this line (inches):
    echo "set terminal pdf size 10, 6" > ${TMP}
    echo "set output \"${DIR}/hist/${FILE}.pdf\"" >> ${TMP}
    if [ $# -eq 2 ]; then
        echo "set key off" >> ${TMP}
    else
    	echo "set key outside" >> ${TMP}
        #echo "set key below samplen 1" >> ${TMP}
    fi
    echo "set title \"${FILE} Histogram\"" >> ${TMP}
    echo "set xlabel \"Depth Complexity\" " >> ${TMP}
    echo "set ylabel \"Frequency (in %)\"" >> ${TMP}
		echo "set logscale y" >> ${TMP}
		echo "set grid" >> ${TMP}
		echo "set xtics 2" >> ${TMP}
    #echo "set ytics 10" >> ${TMP}
    #echo "set yrange [ 0.001 : 100 ]" >> ${TMP}
    #echo "set xrange [ 1 : *]" >> ${TMP}
 
    let COUNT=COUNT+1

    #PLOT
    SUM="$(gawk '{y+=$2}; END {print y}' $FILEPATH)"    
    echo -n "plot \"${FILEPATH}\" using 1:(\$2*100/$SUM) every ::1 w lines title 'step=${COUNT}'" >> ${TMP}
    #echo -n "plot \"${FILEPATH}\" using 1:(log10(\$2+1)+1) every ::1 w lines title 'step=${COUNT}'" >> ${TMP}
 
    echo "" >> ${TMP}
    echo "set output \"trash.pdf\"" >> ${TMP}
    gnuplot ${TMP}
    rm ${TMP}
    rm trash.pdf
done

