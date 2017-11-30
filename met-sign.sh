#!/bin/bash
set -eu

# re-extract DNA methylation, at custom sequence contexts, from BAM file and produce cytosine report
# perform in 4_bismark output sub-directory of wgbs workflow

if [ "$#" -ne 6 ]; then
echo "USAGE: met-sign.sh <SE/PE> <context> <file> <annotation> <sample> <outname>"
echo "EXAMPLE: met-sign.sh SE CHH sample.bam Araport_mRNA.sorted.bed sample mRNA"
exit 1
fi

layout=$1
context=$2
fl=$3
annopath=$4
sample=$5
outname=$6

echo "Extracting CX report from $1 BAM ..."

if [ $layout == "SE" ]; then
bismark_methylation_extractor --comprehensive --multicore 4 --cytosine_report --CX --genome_folder ~/TAIR10/  --report --buffer_size 8G -s ${fl}
fi

if [ $layout == "PE" ]; then
bismark_methylation_extractor --comprehensive --multicore 4 --cytosine_report --CX --genome_folder ~/TAIR10/  --report --buffer_size 8G -p ${fl}
fi

gzip -d *cov.gz
bedfile="${fl::-3}bismark.cov"

echo "reports extracted"
echo "$context from $bedfile"
echo $1 $2 $3 $4 $5 $6

# re-organise report, grep context, and awk to remove C and M
awk '{print $1 "\t" $2 "\t" $2+1 "\t" $6 "\t" $7}' ${fl::-3}CX_report.txt | grep "$context" | awk -F$'\t' ' $1 != "ChrC" && $1 != "ChrM" ' > ${fl::-3}${context}_report.bed

# intersect sub-context info to bismark.cov
intersectBed -wo -a ${fl::-3}${context}_report.bed -b $bedfile | awk 'BEGIN { OFS = "\t" } {print $1, $2, $3, $4, $5, $9}' > ${sample}-${outname}-sub${context}-report.bed

# closest to get info across annotation file and subset to within 1kb
closestBed -D "b" -a ${sample}-${outname}-sub${context}-report.bed -b $annopath | awk -F$'\t' '$NF<1000 && $NF>-1000' > ${sample}-${outname}-sub${context}-report.1k.bed

echo "done... cleaning..."

rm C*txt
rm *_report.txt
rm *bedGraph.gz
rm *M-bias.txt
rm $bedfile
rm ${fl::-3}${context}_report.bed
rm ${sample}-${outname}-sub${context}-report.bed

echo "R"

Rscript ~/scripts/rel_methylation_plots_v2.r ${sample} ${outname} ${context}

echo "DONE"
