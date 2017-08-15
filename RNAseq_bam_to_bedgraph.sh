#!/bin/bash
set -e 
set -u
# Produce coverage data from BAM files for RNAseq or ChIP data in bedgraph format
# Run in directory with sam converted, sorted, indexed  bam file
# Ensure subread indexed genome & chromosome sizes are prepared eg TAIR10/subread_index/tair10.sizes.genome
# samtools faidx tair10.fa | cut -f1,2 tair10.fa.fai > tair10.sizes.genome

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: bam_to_bedgraph.sh <sample name> <unstranded, SE_stranded, PE_stranded>"
echo "EXAMPLE: RNAseq_bam_to_bedgraph.sh col0-r1 unstranded"
exit 1
fi

smp="$1.sorted.bam"
lay=$2

# file for length of all 7 chromosomes
chrc_sizes=chrc_sizes=/home/diep/TAIR10/Arabidopsis_thaliana.TAIR10.22.dna.genome.fa.len

echo "sample = $1"
echo "layout = $2"

echo "Produce ${2} tiled data files from ${smp} ..."

if [[ "$lay" == "unstranded" ]]; then
echo "non-stranded bedgraph"
# non-stranded bedgraph
bedtools genomecov -bga -split -ibam $smp -g $chrc_sizes > ${smp%%sorted.bam*}.bedgraph | tee -a ../${1}*.log

# non-stranded bedgraph with splicing & nt resolution
bedtools genomecov -d -split -ibam $smp -g $chrc_sizes > ${smp%%sorted.bam*}.bed | tee -a ../${1}*.log

fi

if [[ "$lay" == "SE_stranded" ]];then
echo "stranded SE bedgraph"
echo "to be built"

# https://www.biostars.org/p/179035/

# samtools view -f 128 F 4 -b $smp > 

fi

if [[ "$lay" == "PE_stranded" ]];then
echo "stranded PE bedgraph"
echo "Extract properly-paired reads and their mates (ie flags 99/147/163/83) from paired-end BAM files"
# Extract properly-paired reads and their mates (ie flags 99/147/163/83) from paired-end BAM files
# https://gist.github.com/mtw/7175143
# http://seqanswers.com/forums/showthread.php?t=29399
# stranded PE bams

#R1 forward strand
samtools view -f 99 -b ${smp} > ${smp%%bam}R1F.bam 2>&1 | tee -a ../${1}*.log
#R2 reverse strand
samtools view -f 147 -b ${smp} > ${smp%%bam}R2R.bam 2>&1 | tee -a ../${1}*.log
#FORWARD '+' reads
samtools merge -f ${smp%%bam}forward.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam 2>&1 | tee -a ../${1}*.log

#R1 reverse strand
samtools view -f 83 -b ${smp} > ${smp%%bam}R1R.bam | tee -a ../${1}*.log
#R2 forward strand
samtools view -f 163 -b ${smp} > ${smp%%bam}R2F.bam | tee -a ../${1}*.log
# REVERSE '-' reads
samtools merge -f ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam | tee -a ../${1}*.log
# clean
rm ${smp%%bam}*.R*.bam

fi
