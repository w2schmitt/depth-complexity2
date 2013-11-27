#!/bin/bash

DIR=$1
DIR_OUTPUT=$2

for file in $(ls $DIR)
do
	echo `sed -e's/[0-9]*\://g' -e's/{.*}/\[&\]/' -e's/[{}]//g' $DIR/$file > $DIR_OUTPUT/$file`
done 