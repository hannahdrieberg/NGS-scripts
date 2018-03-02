#!/bin/bash
set -u
# Obtaining single mC information from bed files across DMRs identified in `100bp_dmrs.v0.1.sh`
# Handy when throwing lots of samples against DMRs

if [ "$#" -ne 2 ]; then
echo "USAGE: 100bp_heatmap.sh <context> <100bp_dmrs output>"
exit 1
fi

context=$1
dmrs=$2

echo ${context}
echo ${dmrs}

echo "Performing intersectBed..."

for FILE in *_${context}.bed
do
	intersectBed -wo -a ${dmrs} -b $FILE > DMRs-${FILE}
	bedtools groupby -i DMRs-${FILE} -g 1,2,3 -c 9 -o mean > avgDMRs-${FILE}
done
