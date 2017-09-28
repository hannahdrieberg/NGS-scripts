samtools depth  *sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
