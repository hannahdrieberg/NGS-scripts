#!/bin/bash
set -eu

# Produce 100bp windows (100bp.bed) of RNAseq coverage from BAMs across annotations of interest
# Run in directory with sam converted, sorted, indexed  bam file
# Provide path of genome .fa file to produce 100bp genome file

if [ "$#" -lt 4 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_bam_to_100bpwigs.sh <sorted bam> <genome fasta> <annotation> <out>"
echo "EXAMPLE: RNAseq_bam_to_100bpwigs.sh col0-r1 /home/diep/TAIR10/TAIR10_Chr.all.fasta /home/diep/Araport11/annotations/Araport11_TE.bed TE"
exit 1
fi

bam=$1
fas=$2
bedfile=$3
out=$4

echo 'Make 100bp genome bed ...'

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
bedtools coverage -sorted -a temp.genome.100bp.bed -b $bam > ${bam%%.sorted*}_100bp.bed 

echo 'cleaning ...'
# CLEAN
rm temp.genome*

# remove 'Chr'
awk '{print $1=substr($1,4)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ${bam%%.sorted*}_100bp.bed > temp.bed
# sort Bed
sortBed -i temp.bed > ${bam%%.sorted*}_100bp.sorted.bed

echo 'bedtools ...'
# bedtools to desired annotation
closestBed -D "b" -a ${bam%%.sorted*}_100bp.sorted.bed -b $bedfile > ${bam%%.sorted*}_$out.bed

echo 'subset to +1k/-1k ...'
# awk to subset
awk -F$'\t' '$NF<1000 && $NF>-1000' ${bam%%.sorted*}_$out.bed > ${bam%%.sorted*}_$out.1k.bed

rm temp.bed *TE.bed *100bp.bed

