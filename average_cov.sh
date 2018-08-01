#!/bin/bash
set -eu
bam=$1
samtools sort -@ 4 $bam | samtools depth - | awk '{sum+=$3} END { print "Mean depth = ",sum/NR}'
