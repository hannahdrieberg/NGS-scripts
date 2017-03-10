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

scythe -a /home/diep/scripts/TruSeq-adapters.fa -p 0.1 ../$fq1 > ${fq1%%.fastq*}_noadapt.fastq

gzip ${fq1%%.fastq*}_noadapt.fastq

scythe -a /home/diep/scripts/TruSeq-adapters.fa -p 0.1 ../$fq2 > ${fq2%%.fastq*}_noadapt.fastq

gzip ${fq2%%.fastq*}_noadapt.fastq 

sickle pe -g -f ${fq1%%.fastq*}_noadapt.fastq.gz -r ${fq2%%.fastq*}_noadapt.fastq.gz -o ${fq1%%.fastq*}_trimmed.fastq.gz -p ${fq2%%.fastq*}_trimmed.fastq.gz -s trimmed.singles.fastq.gz -t sanger -q 20 -l 20
cd ../

# fastqc again
mkdir 3_trimmed_fastqc
fastqc 2_scythe_sickle/${fq1%%.fastq*}_trimmed.fastq.gz 2_scythe_sickle/${fq2%%.fastq*}_trimmed.fastq.gz

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

# gzip input weirdness ... won't work
# use subread on gzipped paired read files
# subread-align -T 5 -H -t 1 -u -i ${index} --gzFASTQinput -r ${fq1%%.fastq*}_trimmed.fastq.gz -R ${fq2%%.fastq*}_trimmed.fastq.gz -o "${fileID}.sam"

# align un-zipped read pair files
if [[ ${fq1%%.fastq*}_trimmed.fastq.gz == *".gz" ]]; then gzip -d ${fq1%%.fastq*}_trimmed.fastq.gz; fi
if [[ ${fq2%%.fastq*}_trimmed.fastq.gz == *".gz" ]]; then gzip -d ${fq2%%.fastq*}_trimmed.fastq.gz; fi

subread-align -T 5 -H -t 1 -u -i ${index} -r ${fq1%%.fastq*}_trimmed.fastq -R ${fq2%%.fastq*}_trimmed.fastq -o "${fileID}.sam"

# sam to bam
samtools view -S -u $outsam > ${tmpbam}
samtools sort -m 2G -o $outbam ${tmpbam}
samtools index ${outbam}

samtools view -b -S -h ${fq_file%%.fastq*}_trimmed*.sam > ${fq_file%%.fastq*}_trimmed.fq_bismark.bam
samtools sort ${fq_file%%.fastq*}_trimmed.fq_bismark.bam ${fq_file%%.fastq*}_trimmed.fq_bismark.sorted 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools index ${fq_file%%.fastq*}_trimmed.fq_bismark.sorted.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log


