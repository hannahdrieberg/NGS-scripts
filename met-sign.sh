#!/bin/bash
set -e

# re-extract DNA methylation, at custom sequence contexts, from sam file and produce cytosine report
# perform in 4_bismark output sub-directory of wgbs workflow

if [ "$#" -ne 6 ]; then
echo "USAGE: met-sign.sh <context> <file> <file path to met bedfile> <annotation file> <sample> <outname>"
exit 1
fi

context=$1
fl=$2
bedfile=$3
annopath=$4
sample=$5
outname=$6

echo "Extracting CX report from SAM ..."

bismark_methylation_extractor --comprehensive --cytosine_report --CX --genome_folder ~/TAIR10_bs/  --report --buffer_size 8G -s ${fl}

echo "done"
echo "grepping & bedtools ..."

if [ ${context} = CHH ]; then
seq="CAA CAC CAT CCA CCC CCT CTA CTC CTT"
fi

echo ${context}
echo ${fl}
echo ${bedfile}
echo ${annopath}
echo ${sample}
echo ${outname}
echo ${seq}

awk '{print $1 "\t" $2 "\t" $2+1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' ${fl::-3}.CX_report.txt > ${fl::-3}.CX_report.bed

# use grep to get output of specific sequence context from report files

for FILE in $seq
do
grep ${FILE} ${fl::-3}.CX_report.bed > ${FILE}.bed
intersectBed -wo -a ${FILE}.bed -b ${bedfile} > ${context}-${FILE}.bed
sort -k1,1 -k2,2n ${context}-${FILE}.bed -o sorted-${context}-${FILE}.bed
closestBed -D "ref" -a sorted-${context}-${FILE}.bed -b ${annopath} > ${outname}-${context}-${FILE}.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${outname}-${context}-${FILE}.bed > ${outname}-${context}-${FILE}.1k.bed
mv ${outname}-${context}-${FILE}.1k.bed ${sample}-${outname}-${context}-${FILE}.1k.bed
rm ${FILE}.bed
rm ${context}-${FILE}.bed
rm sorted-${context}-${FILE}.bed
rm ${outname}-${context}-${FILE}.bed
done

echo "done"
echo "cleaning ..."

rm *bedGraph.gz
rm *cov.gz
rm *context*bismark.sam.gz.txt
rm *M-bias.txt
rm *splitting_report.txt

echo "done ... r plotting"

Rscript ~/scripts/rel_methylation_plots-v2.r ${sample} ${outname} ${context}

rm *1k.bed
rm *CX*

echo "DONE"
