#!/bin/bash
set -e
set -u

# Use aligned and indexed reads from RNAseq_v0.1.sh to perform DEG analysis 
# based on pedrocrisp/NGS-pipelines/RNAseqPipe3

# feature Counts
# make sure to have genome size file
# samtools faidx tair10.fa
# cut -f1,2 tair10.fa.fai > tair10.sizes.genome


