#!/bin/bash
# Parameters: 
# $1 -> number of polynomials for each degree
# $2 -> m; field size is 2^m
# Degree range from $3 to $4

QUANTITY=$1
M=$2
BEG=$3
END=$4

rm p
./build $M

for i in {$BEG..$END}
do
	echo "$i..."
	./p $QUANTITY $M $i
done