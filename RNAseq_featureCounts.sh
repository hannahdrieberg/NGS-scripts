#!/bin/bash
# Use featureCounts to perform feature counting on BAM files
set -ex

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_featureCounts.sh <filename> <SE/PE> <0 (unstranded)/ 1(stranded) / 2 (reverse stranded)> <bedfile> <outname>"
echo "EXAMPLE: RNAseq_featureCounts.sh col0-r1 PE 1 /home/diep/Araport11/annotations/Araport11_mRNA.bed mRNA"
exit 1
fi

sample="$1.sorted.bam"
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

if [ $layout == "SE" ] ; then 

featureCounts\
	-F SAF\
	-C\
	-T 2\
	-s $strand\
	-a $bedfile\
	-o "${sample}_${outname}.counts"\
	"${sample}"
fi
	
if [ $layout == "PE" ] ; then 

featureCounts\
	-F SAF\
	-p\
	-C\
	-T 2\	
	-s $strand;
	-a $bedfile;
	-o "${sample}_${outname}.counts"\
	"${sample}"
fi

echo "DONE"
