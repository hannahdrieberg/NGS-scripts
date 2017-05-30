#!/bin/bash

# use on sorted BAM file output from bismark alignment to call SNPs using default BS-SNPer settings

## require: 
# perl
# BS-SNPer http://bioinformatics.oxfordjournals.org/content/31/24/4006.long

if [ "$#" -ne 1 ]; then
echo "USAGE: <file>_trimmed.fq_bismark.sorted.bam"
echo "EXAMPLE: BS-SNPer.sh 277_D9"
exit 1
fi

file=$1

perl ~/bin/BS-Snper-master/BS-Snper.pl --fa ~/TAIR10_bs/Arabidopsis_thaliana.TAIR10.22.dna.genome.fa --input ${file}_trimmed.fq_bismark.sorted.bam --output temp.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.15 --minhomfreq 0.85 --minquali 30 --mincover 15 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 > ${file}.SNP.bed 2>${file}_ERR.log

sort -k1,1 -k2,2n ${file}.SNP.bed > ${file}.SNP.bed

