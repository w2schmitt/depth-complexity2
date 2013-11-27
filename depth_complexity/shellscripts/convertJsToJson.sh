#!/bin/bash

dir1=$1
dir2=$2

for file in $(ls $dir1) 
do
	#f="$(basename ${file})"
	name=${file%%.*}.json
	var=$(grep -oE "\-?[0-9]+\.[0-9]+,?]?" $dir1$file | tr -d "\n")
	echo -e "{\n\"frequency\":[$var\n}" > $dir2$name
done
