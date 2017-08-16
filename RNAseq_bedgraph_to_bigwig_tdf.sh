#!/bin/bash
set -e 
set -u
# Use bedgraph files from RNAseq_bam_to_bedgraph.sh to produce bigwigs or TDFs for viewing delight 
# Run in directory with bedgraph filese
# Ensure subread indexed genome & chromosome sizes are prepared eg TAIR10/TAIR10_Chr.all.fasta.len
# samtools faidx TAIR10_Chr.all.fasta | cut -f1,2 TAIR10_Chr.all.fasta.fai > TAIR10_Chr.all.fasta.len

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_bedgraph_to_bigwig_tdf.sh <bedgraph> <unstranded, SE_stranded, PE_stranded> <igv genome>"
echo "EXAMPLE: RNAseq_bedgraph_to_bigwif_tdf.sh col0-r1 unstranded /home/diep/araport11_igv_genome/Araport11.genome"
exit 1
fi

bg=$1
lay=$2
igv_genome=$3

# file for length of all 7 chromosomes
chrc_sizes=chrc_sizes=/home/diep/TAIR10/TAIR10_Chr.all.fasta.len

echo "sample = $1"
echo "layout = $2"

echo "Produce ${2} tiled data files from ${smp} ..."

if [[ "$lay" == "unstranded" ]]; then
echo "non-stranded bedgraph"
# non-stranded bedgraph
bedtools genomecov -bga -split -ibam $smp -g $chrc_sizes > ${smp%%sorted.bam*}.bedgraph
fi

if [[ "$lay" == "SE_stranded" ]];then
echo "blah"
fi

if [[ "$lay" == "PE_stranded" ]];then
echo "Extract properly-paired reads and their mates (ie flags 99/147/163/83) from paired-end BAM files"

fi
