#!/bin/bash

# Obtaining single mC information from bed files across DMRs identified in `100bp_dmrs.v0.1.sh`
# Handy when throwing lots of samples against DMRs

echo "USAGE: 100bp_heatmap.sh <context> <100bp_dmrs output>"

context=$1
dmrs=$2

echo "Performing intersectBed..."

for FILE in *_{context}.bed
do
	intersectBed -wo -a ${dmrs} -b $FILE > DMRs-${FILE}.bed
	bedtools groupby -i DMRs-${FILE} -g 1,2,3 -c 9 -o mean > DMRs-${FILE}-sum.bed
done