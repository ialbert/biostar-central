#!/bin/bash

YEARS=`seq 2001 2014`
MONTHS=(January February March April May June July August September October November December)

COUNTER=1000

for y in ${YEARS[@]}; do
	for m in ${MONTHS[@]}; do
		url="https://stat.ethz.ch/pipermail/bioconductor/$y-$m.txt.gz"
		out="$COUNTER-$y-$m.txt.gz"
		wget $url -O $out
		COUNTER=$((COUNTER + 1))
	done
done