#!/bin/bash
# example: bash exec_polys.sh -f polynomials.txt -e main -c 0

while getopts e:f:c: flag
do
    case "${flag}" in
        f) filename=${OPTARG};;
		e) execut=${OPTARG};;
		c) count=${OPTARG};;
    esac
done;

M=$(head -n 1 $filename)
T=$(sed -n '2p' < $filename)

{
	for ((i=2+($count);i--;));
	do read; done;

	EXPECT=$T
	# EXPECT=128

	while read line;
	do
		IFS=' ' read -a params <<< "$line"
		LINE=( $line )
		length=${#LINE[@]}
		DIFF=$((EXPECT+1-length))
		zeros=$(for ((i=$DIFF;i--;)) ; do echo "0"; done)

	    params=$(echo $line $zeros $M $EXPECT)
	    ./$execut $params
	done;
	
} < $filename