#!/bin/bash
set -eu

# Produce 100bp windows (100bp.bed) of RNAseq coverage from BAMs
# Run in directory with sam converted, sorted, indexed  bam file
# Provide path of genome .fa file to produce 100bp genome file

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_bam_to_100bpwigs.sh <bam file> <genome fasta>"
echo "EXAMPLE: RNAseq_bam_to_100bpwigs.sh col0-r1 /home/diep/TAIR10/TAIR10_Chr.all.fasta"
exit 1
fi

bam=$1
fas=$2

echo 'Make 100bp genome bed'

# use samtools to generate fasta index
samtools faidx $fas
# use awk on index to make genome file
# https://www.biostars.org/p/70795/
awk -v OFS='\t' {'print $1,$2'} ${fas}.fai > temp.genome 
# use genome file to make 100bp windows across genome
bedtools makewindows -g temp.genome -w 100 > temp.genome.100bp.bed
sort -k1,1 -k2,2n temp.genome.100bp.bed > temp.genome.100bp.sorted.bed

# use bedtools coverage to get coverage across 100bp windows from BAM
# MAKE SURE TO USE -sorted FLAG
bedtools coverage -sorted -a temp.genome.100bp.bed -b $bam > ${bam%%sorted}100bp.bed 

echo 'cleaning ...'
# CLEAN
rm temp.genome*bed

