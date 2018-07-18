#!/bin/bash

# Assign (feature)Counts to exons based on the AtRTD2 reference transcript dataset from aligned BAM files
# Additional info https://www.biostars.org/p/321379/
# Calixto, C.P.G., Guo, W., James, A.B., Tzioutziou, N.A., Entizne, J.C., Panter, P.E., Knight, H., Nimmo, H., Zhang, R., and Brown, J.W.S. (2018). Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell: tpc.00177.2018.
# https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf 

# To do: include use of either RTD2 or RTD2_quasi

set -eu

if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_featureCounts.sh <filename> <SE/PE> <0/1/2>"
echo "EXAMPLE: RNAseq_featureCounts.sh col0-r1.sorted.bam PE 1"
echo "0 = unstranded, 1 = stranded, 2 = reverse stranded"
exit 1
fi

sample=$1
layout=$2
strand=$3
bedfile="$HOME/AtRTD2/AtRTD2_19April2016.gtf"

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo "$bedfile"
echo ""
echo "Alternative splicing - $layout $strand featureCounts on $bedfile in $sample ..."
echo ""

if [[ $layout == "SE" ]]; then 

featureCounts\
	-F GTF\
	-C\
	-T 4\
	-f\ 	# Perform read counting at feature level (eg. exons rather than genes).
	-t exon\
	-g gene_id\
	-O\     # or you will lose all the reads or read-pairs that overlap with multiple exons
	-M\	# or all the reads that have 2 or more mapping locations are not counted at all
	-s $strand\
	-a $bedfile\
	-o "${1%%.bam*}_RTD2.counts"\
	$sample 2>&1 | tee -a ../*log
fi
	
if [[ $layout == "PE" ]]; then 

featureCounts\
	-F GTF\
	-p\
	-C\
	-T 4\
	-f\
	-t exon\
	-g gene_id\
	-O\
	-M\
	-s $strand\
	-a $bedfile\
	-o "${1%%.bam*}_RTD2.counts"\
	$sample 2>&1 | tee -a ../*log
fi

echo "DONE"
