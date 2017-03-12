#!/bin/bash

# Take raw ChIP-seq reads, align and produce files to view distribution of histone marks in IGV
# Will eventually introduce peak caller
# H3K9me2 for now
# based on pedrocrisp/NGS-pipelines/RNAseqPipe3

# exit if error and echo steps
set -e
set -x

# make sure genome index has been built
# build index
# subread-buildindex -o TAIR10_subread_index TAIR10_chr1.fas TAIR10_chr2.fas TAIR10_chr3.fas TAIR10_chr4.fas TAIR10_chr5.fas TAIR10_chrC.fas TAIR10_chrM.fas

if [ "$#" -lt 5 ]; then
echo "Missing required arguments!"
echo "USAGE: chip-seq_v0.1.sh <SE, PE> <fastq R1> <R2> <ref genome> <indexed genome> <fileID output>"
exit 1
fi

###
# SINGLE END
###

if [ "$1" == "SE" ] then

# require
if [ "$#" -ne 4 ]; then
echo "Missing required arguments for single-end!"
echo "USAGE: chip-seq_v0.1.sh <-se> <R1> <ref genome> <indexed ref> <fileID output>"
exit 1
fi

fi

#### 
## PAIRED END
####

if [ "$1" == "PE" ] then

if [ "$#" -ne 6 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: chip-seq_v0.1.sh <PE> <R1> <R2> <ref genome> <indexed genome> <fileID output>"
exit 1
fi

#gather input variables
type=$1
fq1=$2;
fq2=$3;
genome=$4; #path to genome
index=$5; #path to genome index
fileID=$6;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing ChIP-seq alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq1 $fq2"
echo "genome: $genome"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_ChIP_${dow}
mv $fq1 ${fileID}_ChIP_${dow}
mv $fq2 ${fileID}_ChIP_${dow}
cd ${fileID}_ChIP_${dow}

if [[ $fq1 != *.gz ]];then
        gzip $fq1 > $fq1
fi

if [[ $fq2 != *.gz ]];then
        gzip $fq2 > $fq2
fi

# initial fastqc
mkdir 1_fastqc
fastqc $fq1 $fq2
mv ${fq1%%.fastq*}_fastqc* 1_fastqc
mv ${fq2%%.fastq*}_fastqc* 1_fastqc

# adapter and quality trimming using scythe and sickle 
mkdir 2_scythe_sickle
cd 2_scythe_sickle

scythe -a /home/diep/scripts/TruSeq-adapters.fa -p 0.1 ../$fq1 > ${fq1%%.fastq*}_noadapt.fastq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

scythe -a /home/diep/scripts/TruSeq-adapters.fa -p 0.1 ../$fq2 > ${fq2%%.fastq*}_noadapt.fastq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

sickle pe -g -f ${fq1%%.fastq*}_noadapt.fastq.gz -r ${fq2%%.fastq*}_noadapt.fastq -o ${fq1%%.fastq*}_trimmed.fastq -p ${fq2%%.fastq*}_trimmed.fastq -s trimmed.singles.fastq -t sanger -q 20 -l 20 2>&1 | tee -a ../${fileID}_logs_${dow}.log

rm ${fq1%%.fastq*}_noadapt.fastq
rm ${fq2%%.fastq*}_noadapt.fastq
cd ../

# fastqc again
mkdir 3_trimmed_fastqc
fastqc 2_scythe_sickle/${fq1%%.fastq*}_trimmed.fastq 2_scythe_sickle/${fq2%%.fastq*}_trimmed.fastq

mv 2_scythe_sickle/${fq1%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc
mv 2_scythe_sickle/${fq2%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc

mkdir 0_fastq
mv $fq1 0_fastq
mv $fq2 0_fastq

# subread align
mkdir 4_subread-align
mv 2_scythe_sickle/${fq1%%.fastq*}_trimmed.fastq.gz -t 4_subread-align/
mv 2_scythe_sickle/${fq2%%.fastq*}_trimmed.fastq.gz -t 4_subread-align/
cd 4_subread-align/

echo "Beginning alignment ..."

subread-align -T 5 -H -t 1 -u -i ${index} -r ${fq1%%.fastq*}_trimmed.fastq -R ${fq2%%.fastq*}_trimmed.fastq -o "${fileID}.sam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

if [[ $fq1%%.fastq}* != *".gz" ]]; then gzip ${fq1%%.fastq*}_trimmed.fastq; fi
if [[ $fq2%%.fastq}* != *".gz" ]]; then gzip ${fq2%%.fastq*}_trimmed.fastq; fi

echo "Alignment complete ... making sorted bam file with index ..."

# samtools view to convert the sam file to bam file
tmpbam="${fileID}.temp.bam"
outbam="${fileID}.sorted"

samtools view -S -u ${fileID}.sam > ${tmpbam}

# Sort the temporary bam file by chromosomal position, and save the sorted file.
samtools sort -m 2G ${tmpbam} $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
# Make an index of the sorted bam file
samtools index ${outbam}.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# delete temp bam and gzip sam
rm -v ${tmpbam}
gzip ${fileID}.sam
mv *trimmed.fastq.gz ../2_scythe_sickle/

# feature Counts
# make sure to have genome size file 
# samtools faidx tair10.fa
# cut -f1,2 tair10.fa.fai > tair10.sizes.genome

# could probably make separate file e.g. bam_to_TDF.sh
echo "Produce tiled data files from bam ..."

# bam to tdf
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

echo "Bedgraph"

# tdf from stranded bedgraphs
mkdir stranded-tdf/
mv ${outbam}.forward.bam ${outbam}.reverse.bam -t stranded-tdf/
cd stranded-tdf/

# plus strand
bedtools genomecov -bga -split -ibam ${outbam}.reverse.bam -g /home/diep/TAIR10/subread_index/tair10.sizes.genome > ${outbam}.plus.bedgraph

# minus strand
bedtools genomecov -bga -split -ibam ${outbam}.forward.bam -g /home/diep/TAIR10/subread_index/tair10.sizes.genome > ${outbam}.minus.bedgraph

# make tdf
igvtools toTDF ${outbam}.plus.bedgraph ${outbam}.plus.tdf /home/diep/araport11_igv_genome/Araport11.genome
igvtools toTDF ${outbam}.minus.bedgraph ${outbam}.minus.tdf /home/diep/araport11_igv_genome/Araport11.genome



