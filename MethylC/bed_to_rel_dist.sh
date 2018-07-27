#!/bin/bash
# Plotting DNA methylation over gene models

#######################
# REQUIREMENTS
# bedtools
# awk
# R with libraries: fields
#######################

if [ "$#" -ne 3 ]; then
echo "USAGE: bed_to_rel_dist.sh <input path to bed file> <filename prefix> <output map name>"
echo "EXAMPLE: bed_to_rel_dist.sh $HOME/Araport11/Araport11_genes.sorted.bed sample-r1 genes"
exit 1
fi

bedpath=$1
filename=$2
outname=$3

sort -k1,1 -k2,2n ${filename}_CG_100bp*.bed -o ${filename}_CG_100bp.bed
sort -k1,1 -k2,2n ${filename}_CHG_100bp*.bed -o ${filename}_CHG_100bp.bed
sort -k1,1 -k2,2n ${filename}_CHH_100bp*.bed -o ${filename}_CHH_100bp.bed

#get total number of columns for both input files
l1="$(cat ${filename}_CG_100bp.bed | awk 'BEGIN{FS="\t"};{print NF}' | head -n 1)"
l2="$(cat ${bedpath} | awk 'BEGIN{FS="\t"};{print NF}' | head -n 1)"

#convert the wigs to bed
####################### some bedtools stuff
echo "Performing closestBed of CHG methylation..."
closestBed -D "ref" -a ${filename}_CHG_100bp.bed -b $bedpath > ${filename}_CHG_${outname}.bed
echo "Performing closestBed of CHH methylation..."
closestBed -D "ref" -a ${filename}_CHH_100bp.bed -b $bedpath > ${filename}_CHH_${outname}.bed
echo "Performing closestBed of CG methylation..."
closestBed -D "ref" -a ${filename}_CG_100bp.bed -b $bedpath > ${filename}_CG_${outname}.bed

#subset to the regions within 100bp of a gene (make the files more manageable for R)
echo "subsetting files to within 1kb..."
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CHG_${outname}.bed > ${filename}_CHG_${outname}.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CHH_${outname}.bed > ${filename}_CHH_${outname}.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CG_${outname}.bed > ${filename}_CG_${outname}.1k.bed
#######################

rm ${filename}_CG_100bp.bed ${filename}_CHG_100bp.bed ${filename}_CHH_100bp.bed
rm ${filename}_CHG_${outname}.bed ${filename}_CHH_${outname}.bed ${filename}_CG_${outname}.bed

echo "Performing R plots..."
#initiate the R script to create the plots
Rscript $HOME/scripts/MethylC/rel_methylation_plots.r ${filename} ${outname} ${l1} ${l2}
