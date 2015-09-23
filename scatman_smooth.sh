#!/bin/bash

# Observe whether annotation feature length correlates with methylation level (technical bias).
# E.g. Do longer TEs have higher levels of methylation? 
# intersect bed file with annotation file, then produce scattersmooth plot of methylation level by length of feature.

if [ "$#" -ne 3 ]; then
echo "scatman_smooth.sh annotation sample outname"
echo "e.g. scatman_smooth.sh ./TAIR_TE_subset.bed 317-1-4 TE" 
exit 1
fi

annotation=$1	
sample=$2
outname=$3

# intersect bed files with annotation file
echo "Performing intersectBed of CHG methylation..."
intersectBed -wa -wb -a ${filename}_CHG_100bp.bed -b $bedpath > ${filename}_CHG_${outname}.bed
echo "Performing intersectBed of CHH methylation..."
intersectBed -wa -wb -a ${filename}_CHH_100bp.bed -b $bedpath > ${filename}_CHH_${outname}.bed
echo "Performing intersectBed of CpG methylation..."
intersectBed -wa -wb -a ${filename}_CpG_100bp.bed -b $bedpath > ${filename}_CpG_${outname}.bed

# take output and produce smoothed scatterplots
echo "Performing scatterSmooth in R"

Rscript $HOME/scripts/smooth_scat.r ${filename} ${outname}