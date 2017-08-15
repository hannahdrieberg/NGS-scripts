#!/bin/bash

# Produce tiled data files for RNAseq or ChIP data from sorted BAMs for viewing delight
# Run in directory with sam converted, sorted, indexed  bam file

set -u  

if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: bam_to_tdf.sh <sample name> <unstranded, SE_stranded, PE_stranded> <igv genome path>"
echo "eg. bam_to_tdf.sh col0-H3K9me2-r1 unstranded /home/diep/araport11_igv_genome/Araport11.genome"
exit 1
fi

smp="$1.sorted.bam"
lay=$2
igvnome=$3

echo "sample = $1"
echo "layout = $2"
echo "igv genome path = $3"

echo "Produce ${2} tiled data files from ${smp} ..."

if [[ "$lay" == "unstranded" ]]; then

echo ""
echo "tdf from non-stranded bedgraph"
echo ""

# non-stranded bedgraph
bedtools genomecov -bga -split -ibam $smp -g /home/diep/TAIR10/subread_index/tair10.sizes.genome > ${smp%%sorted.bam*}.bedgraph

# make tdf
igvtools toTDF ${smp%%sorted.bam*}.bedgraph ${1}.tdf ${igvnome}

rm ${smp%%sorted.bam*}.bedgraph

fi

if [[ "$lay" == "SE_stranded" ]];then
echo "blah"
fi

if [[ "$lay" == "PE_stranded" ]];then

echo "Extract properly-paired reads and their mates (ie flags 99/147/163/83) from paired-end BAM files"

# Extract properly-paired reads and their mates (ie flags 99/147/163/83) from paired-end BAM files
# https://gist.github.com/mtw/7175143
# http://seqanswers.com/forums/showthread.php?t=29399
# make sure there are indexed chromosome files with samtools faidx

samtools view -b -f99 "${outbam}.bam" > ${outbam}.R1F.bam
samtools view -b -f147 "${outbam}.bam" > ${outbam}.R2R.bam
samtools merge -f ${outbam}.forward.bam ${outbam}.R1F.bam ${outbam}.R2R.bam

samtools view -b -f163 "${outbam}.bam" > ${outbam}.R2F.bam
samtools view -b -f83 "${outbam}.bam" > ${outbam}.R1R.bam
samtools merge -f ${outbam}.reverse.bam ${outbam}.R1R.bam ${outbam}.R2F.bam

rm ${outbam}*.R*.bam

echo ""
echo "tdf from stranded bedgraphs"
echo ""

# tdf from stranded bedgraphs
mkdir stranded-tdf/
mv ${outbam}.forward.bam ${outbam}.reverse.bam -t stranded-tdf/
cd stranded-tdf/

# plus strand
bedtools genomecov -bga -split -ibam ${outbam}.reverse.bam -g /home/diep/TAIR10/subread_index/tair10.sizes.genome > ${outbam}.plus.bedgraph

# minus strand
bedtools genomecov -bga -split -ibam ${outbam}.forward.bam -g /home/diep/TAIR10/subread_index/tair10.sizes.genome > ${outbam}.minus.bedgraph

# make tdf
igvtools toTDF ${outbam}.plus.bedgraph ${outbam}.plus.tdf ${igvnome}
igvtools toTDF ${outbam}.minus.bedgraph ${outbam}.minus.tdf ${igvnome}

rm ${outbam}.plus.bedgraph
rm ${outbam}.minus.bedgraph

fi


