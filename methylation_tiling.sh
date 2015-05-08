#!/bin/bash

# Script similar to bed_to_rel, however, in this case producing a file with methylation values at features
# of interest identified with an annotation file. e.g. Methylation near and across all protein coding genes.
# Provide the script WIG files of interest and an annotation file.

if [ "$#" -ne 4 ]; then
##############################################################################################################
echo "USAGE: methylation_tiling.sh <filename prefix> <relative path to annotation file> <output map name> <subset>"
echo "Example: methylation_tiling.sh 317-1-4 ../annotation/TAIR10_TE.bed TE 3000"
###############################################################################################################
exit 1
fi

sample=$1
annotation=$2
outname=$3
subset=$4

sed -e "1d" ${sample}_CpG_100bp.wig > ${sample}_CpG_100bp.bed
sed -e "1d" ${sample}_CHG_100bp.wig > ${sample}_CHG_100bp.bed
sed -e "1d" ${sample}_CHH_100bp.wig > ${sample}_CHH_100bp.bed

echo "Performing closestBed of CHG methylation..."
closestBed -D "ref" -a ${sample}_CHG_100bp.bed -b ${annotation} > ${sample}_CHG_${outname}.bed
echo "Performing closestBed of CHH methylation..."
closestBed -D "ref" -a ${sample}_CHH_100bp.bed -b ${annotation} > ${sample}_CHH_${outname}.bed
echo "Performing closestBed of CpG methylation..."
closestBed -D "ref" -a ${sample}_CpG_100bp.bed -b ${annotation} > ${sample}_CpG_${outname}.bed

echo "subsetting files"
awk -v var=${subset} -F$'\t' '$NF<var && $NF>-var' ${sample}_CHG_${outname}.bed > ${sample}_CHG_${outname}.${subset}.bed
awk -v var=${subset} -F$'\t' '$NF<var && $NF>-var' ${sample}_CHH_${outname}.bed > ${sample}_CHH_${outname}.${subset}.bed
awk -v var=${subset} -F$'\t' '$NF<var && $NF>-var' ${sample}_CpG_${outname}.bed > ${sample}_CpG_${outname}.${subset}.bed

# Create directory and cleanup

mkdir ${sample}_${outname}
mv ${sample}_CpG_${outname}.${subset}.bed ${sample}.${outname}
mv ${sample}_CHH_${outname}.${subset}.bed ${sample}.${outname}
mv ${sample}_CHG_${outname}.${subset}.bed ${sample}.${outname}

rm ${sample}*.bed
