#!/bin/bash

M=13

rm main
./build $M

for i in {65..95}
do
echo "$i ..."
sh exec_polys.sh -f polynomials/polys-2500-$M-$i.txt -e main -c 0 >> results-2500-$M-$i.txt
done