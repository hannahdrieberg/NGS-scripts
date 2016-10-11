#!/bin/bash

perl ~/bin/BS-Snper-master/BS-Snper.pl --fa ~/TAIR10_bs/Arabidopsis_thaliana.TAIR10.22.dna.genome.fa --input sample.fq_bismark.sorted.bam --output temp.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.15 --minhomfreq 0.85 --minquali 30 --mincover 15 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 > sample_SNP.out 2>ERR.log
