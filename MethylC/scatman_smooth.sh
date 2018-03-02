#!/bin/bash

# Observe whether annotation feature length correlates with methylation level (technical bias).
# E.g. Do longer TEs have higher levels of methylation? 
# Intersect bed file with annotation file;
# Produce scattersmooth plot of methylation level vs feature length.

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
intersectBed -wa -wb -a ${sample}_CHG.bed -b $annotation > ${sample}_CHG_${outname}.bed
echo "Performing intersectBed of CHH methylation..."
intersectBed -wa -wb -a ${sample}_CHH.bed -b $annotation > ${sample}_CHH_${outname}.bed
echo "Performing intersectBed of CpG methylation..."
intersectBed -wa -wb -a ${sample}_CpG.bed -b $annotation > ${sample}_CpG_${outname}.bed

# take output and produce smoothed scatterplots
echo "Performing scatterSmooth in R"

Rscript $HOME/scripts/smooth_scat.r ${sample} ${outname}

echo "cleanup intermediates"
rm ${sample}_*_${outname}.bed
