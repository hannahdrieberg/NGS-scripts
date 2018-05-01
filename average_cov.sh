#!/bin/bash
set -eu

bam=$1
 
samtools depth  $bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
