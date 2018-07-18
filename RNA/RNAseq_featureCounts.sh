#!/bin/bash
# Use featureCounts to assign counts to annotated features from aligned BAM files

set -eu

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_featureCounts.sh <filename> <SE/PE> <0/1/2> <bedfile> <outname>"
echo "EXAMPLE: RNAseq_featureCounts.sh col0-r1.sorted.bam PE 1 /home/diep/Araport11/annotations/Araport11_mRNA.bed mRNA"
echo "0 = unstranded, 1 = stranded, 2 = reverse stranded"
exit 1
fi

sample=$1
layout=$2
strand=$3
bedfile=$4
outname=$5

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "bedfile = $4 ($5)"
echo ""
echo "$layout $strand featureCounts on $bedfile ($outname) in $sample ..."
echo ""

# SAF formatted file from BED file anno
awk -F'\t' '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $bedfile > temp.saf
awk 'BEGIN {print "GeneID""\t""Chr""\t""Start""\t""End""\t""Strand"}{print}' temp.saf > temp2.saf

if [[ $layout == "SE" ]]; then 

featureCounts\
	-F 'SAF'\
	-C\
	-T 2\
	-s $strand\
	-a temp2.saf\
	-o "${1%%.bam*}_${outname}.counts"\
	$sample 2>&1 | tee -a ../*log
fi
	
if [[ $layout == "PE" ]]; then 

featureCounts\
	-F SAF\
	-p\
	-C\
	-T 2\
	-s $strand\
	-a temp2.saf\
	-o "${1%%.bam*}_${outname}.counts"\
	$sample 2>&1 | tee -a ../*log
fi

rm temp*.saf -v

echo "DONE"
