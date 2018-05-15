#!/bin/bash
set -e 
set -u

# Produce coverage data from BAM files for RNAseq or ChIP data in bedgraph format
# Then produce bigWigs files for viewing delight
# Run in directory with sam converted, sorted, indexed  bam file
# Ensure subread indexed genome & chromosome sizes are prepared eg TAIR10/TAIR10_Chr.all.fasta.len
# samtools faidx TAIR10_Chr.all.fasta | cut -f1,2 TAIR10_Chr.all.fasta.fai > TAIR10_Chr.all.fasta.len

if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: RNAseq_bam_to_bedgraph.sh <sample name> <SE, PE> <unstranded, stranded>"
echo "EXAMPLE: RNAseq_bam_to_bedgraph.sh col0-r1.bam PE unstranded"
exit 1
fi

smp=$1
lay=$2
str=$3

# file for length of all 7 chromosomes
chrc_sizes=/home/diep/TAIR10/TAIR10_Chr.all.fasta.len

echo ""
echo "sample = $1"
echo "layout = $2"
echo "strand = $3"
echo ""
echo "Produce $lay $str bigWig file(s) from $smp ..."
echo ""

if [[ "$lay" == "SE" ]] && [[ "$str"  == "unstranded" ]] ; then 

echo "BAM to bedgraph ..."
# unstranded bedgraph
bedtools genomecov -bga -split -ibam $smp -g $chrc_sizes > ${smp%%sorted.bam*}bg

# bg to bigWig
echo "bigWig ..."
/home/diep/bin/kentUtils/bin/bedGraphToBigWig ${smp%%sorted.bam*}bg ${chrc_sizes} ${smp%%bam}bigWig

fi


if [[ "$lay" == "SE" ]] && [[ "$str"  == "stranded" ]] ; then

echo "Does not exist ....bye"
# https://www.biostars.org/p/179035/
# samtools view -f 128 F 4 -b $smp > 

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "unstranded" ]] ; then

# need sorted bam
echo "sort by position"
samtools sort ${smp} -o ${smp%%bam}sorted.bam

echo "BAM to bedgraph ..."
# unstranded bedgraph
bedtools genomecov -bga -split -ibam $smp -g $chrc_sizes > ${smp%%sorted.bam*}bg

# bg to bigWig
echo "bigWig ..."
/home/diep/bin/kentUtils/bin/bedGraphToBigWig ${smp%%sorted.bam*}bg ${chrc_sizes} ${smp%%bam}bigWig

fi


if [[ "$lay" == "PE" ]] && [[ "$str"  == "stranded" ]] ; then

echo "Extract properly-paired read mates (+ flags 99/147; - flags 83/163) from paired-end BAM files"
# http://seqanswers.com/forums/showthread.php?t=29399

# need sorted bam
echo "sort by position"
samtools sort ${smp} -o ${smp%%bam}sorted.bam

echo "forward strand"
# R1 forward
samtools view -f 99 -b ${smp%%bam}sorted.bam > ${smp%%bam}R1F.bam
# R2 reverse
samtools view -f 147 -b ${smp%%bam}sorted.bam > ${smp%%bam}R2R.bam
# FORWARD R1 read pairs
samtools merge -f ${smp%%bam}forward.bam ${smp%%bam}R1F.bam ${smp%%bam}R2R.bam

echo "reverse strand"
# R1 reverse
samtools view -f 83 -b ${smp} > ${smp%%bam}R1R.bam
# R2 forward
samtools view -f 163 -b ${smp} > ${smp%%bam}R2F.bam
# REVERSE R1 read pairs
samtools merge -f ${smp%%bam}reverse.bam ${smp%%bam}R1R.bam ${smp%%bam}R2F.bam

echo "BAM to stranded bedgraph ..."
# plus strand
bedtools genomecov -bga -split -ibam ${smp%%bam}reverse.bam -g $chrc_sizes > ${smp%%bam}plus.bg
# minus strand
bedtools genomecov -bga -split -scale -1 -ibam ${smp%%bam}forward.bam -g $chrc_sizes > ${smp%%bam}minus.bg

echo "bigWigs..."
$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}plus.bg ${chrc_sizes}  ${smp%%bam}plus.bigWig
$HOME/bin/kentUtils/bin/bedGraphToBigWig ${smp%%bam}minus.bg ${chrc_sizes} ${smp%%bam}minus.bigWig

rm ${smp%%.bam}*R1*bam -v
rm ${smp%%.bam}*R2*bam -v
rm ${smp%%.bam}*forward*bam -v
rm ${smp%%.bam}*reverse*bam -v

fi

# clean up tmps
rm ${smp%%.sorted.bam}*bg -v
