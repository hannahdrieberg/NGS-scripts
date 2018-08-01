#!/bin/bash
set -eu

bam=$1

samtools --threads 4 sort $bam | samtools depth - | awk '{sum+=$3} END { print "Mean depth = ",sum/NR}'
