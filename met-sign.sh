#!/bin/bash
set -e

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
bismark_methylation_extractor --comprehensive --multicore 2 --cytosine_report --CX --genome_folder ~/TAIR10/  --report --buffer_size 8G -s ${fl}
fi

if [ $layout == "PE" ]; then
bismark_methylation_extractor --comprehensive --multicore 2 --cytosine_report --CX --genome_folder ~/TAIR10/  --report --buffer_size 8G -p ${fl}
fi

echo "done"
echo "grepping & bedtools ..."

if [ ${context} = CHH ]; then
seq="CAA CAC CAT CCA CCC CCT CTA CTC CTT"
fi

echo $1 $2 $3 $4 $5 $6

gzip -d *cov.gz

awk '{print $1 "\t" $2 "\t" $2+1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' ${fl::-3}CX_report.txt > ${fl::-3}CX_report.bed

# use grep to extract specific sequence context from report files
echo "grep $seq contexts from $context bedfile"

bedfile="${fl::-3}bismark.cov"
sortBed -i $bedfile > ${bedfile}.sorted.bed

for FILE in $seq
do
grep ${FILE} ${fl::-3}CX_report.bed > ${FILE}.bed
sortBed -i {FILE}.bed > ${FILE}.sorted.bed
intersectBed -wo -a ${FILE}.sorted.bed -b ${bedfile}.sorted.bed > ${context}-${FILE}.bed
closestBed -D "b" -a ${context}-${FILE}.bed -b $annopath > ${outname}-${context}-${FILE}.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${outname}-${context}-${FILE}.bed > ${sample}-${outname}-${context}-${FILE}.1k.bed
rm ${FILE}*.bed
rm ${context}-${FILE}.bed
rm sorted-${context}-${FILE}.bed
done

echo "done... cleaning..."

rm C*txt
rm *splitting_report.txt
rm *bedGraph.gz
rm *M-bias.txt
rm ${bedfile}
rm ${bedfile}.sorted.bed

echo "R"

Rscript ~/scripts/rel_methylation_plots-v2.r ${sample} ${outname} ${context}

rm *1k.bed
rm *CX*

echo "DONE"
