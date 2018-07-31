#!/bin/bash
set -eu

# Generate mean methylation levels into custom bins with desired read depth/coverage from per-site BED files

if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: wgbs_custom_bins.sh <sample> <genome fasta> <coverage> <bin size>"
echo "EXAMPLE: wgbs_custom_bins.sh col0-r1 /home/diep/TAIR10/TAIR10_Chr.all.fasta 15 100"
exit 1
fi

bed=$1
fas=$2
cov=$3
bin=$4
window=$(expr $bin - 1)

echo "Weighted methylation in $bed across $bin bp windows with depth >= $cov ..."

# use samtools to generate fasta index
samtools faidx $fas

# use awk on index to make genome file
# https://www.biostars.org/p/70795/
awk -v OFS='\t' {'print $1,$2'} ${fas}.fai > temp.genome

# use genome file to make 100bp windows across genome
bedtools makewindows -g temp.genome -w ${window} -s ${bin} | sortBed | awk -F$'\t' ' $1 != "ChrC" && $1 != "ChrM" ' > temp.genome.${bin}bp.sorted.bed

# use bedtool intersect and groupBy to get mean methylation levels per bin based on per-site methylation
echo 'Bedtools CG ...'
sort -k1,1 -k2,2n ${bed}_CG.bed.bismark.cov | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 7,8,9 -o mean,sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($5 / ($5+$6)*100 ),$5 = ($5 + $6)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CG_${bin}bp_${cov}cov.bed

echo 'Bedtools CHG ...'
sort -k1,1 -k2,2n ${bed}_CHG.bed.bismark.cov | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 7,8,9 -o mean,sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($5 / ($5+$6)*100 ), $5 = ($5 + $6)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CHG_${bin}bp_${cov}cov.bed

echo 'Bedtools CHH ...'
sort -k1,1 -k2,2n ${bed}_CHH.bed.bismark.cov | bedtools intersect -sorted -wo -a temp.genome.${bin}bp.sorted.bed -b "stdin" | groupBy -g 1,2,3 -c 7,8,9 -o mean,sum,sum | awk -v OFS='\t' '{print $1,$2,$3,$4 = ($5 / ($5+$6)*100 ), $5 = ($5 + $6)}' | awk '{ if ($5 >= '$cov') { print } }' > ${bed}_CHH_${bin}bp_${cov}cov.bed

echo 'cleaning ...'
# CLEAN
rm temp.genome*

