#!/bin/bash

# Take raw ChIP-seq reads, align and produce files to view distribution of histone marks in IGV
# Will eventually introduce peak caller
# H3K9me2 for now
# based on pedrocrisp/NGS-pipelines/RNAseqPipe3

# exit if error and echo steps
set -e
set -x

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: wgbs_se_pipelinev0.4.sh <-pe, -se, -se_epi, or -pese> <in fastq R1> <in fastq R2 (if PE)> <path to bismark genome folder> <fileID for output files>"
exit 1
fi




 
