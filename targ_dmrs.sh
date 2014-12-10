#!/bin/bash
set -u #Treat unset variables as an error when substituting.

# targ_dmrs.V1
# The aim of this script is to take .cov output files after the wgbs_pipeline, which
# contain methylation information still at an individual cytosine resolution, and 
# compare this information to interested genomic features such as genes. To
# do this this script will take the .cov file of interest and use bedtools intersect 
# to compare that to an annotation file, which will contain position information about 
# your genes of interest.

# Execute in directory containing required files, supplying both your file of interest
# and your annotation file.

# How to
if [ "$#" -ne 3 ]; then
echo "USAGE: targ_dmrs.sh <Sample> <Context> <Annotation>"
echo "EXAMPLE: targ_dmrs.sh Pete4 CpG mRNA"
echo "Bed file of interest here: 'Pete4_CpG.bed.bismark.cov'"
echo "Context: CpG, CHH or CHG"
echo "Annotation options: mRNA, TE, miRNA"
exit 1
fi

# Define arguments
sample=$1
context=$2
annotation=$3

# Use Bedtools to intersect the bed file with the annotation file
# bedtools intersect [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

bedtools intersect -wa -wb -a ${sample}_${context}.bed.bismark.cov -b TAIR10_${annotation}_subset.bed > ${sample}_${context}_${annotation}_intersect_output.bed

