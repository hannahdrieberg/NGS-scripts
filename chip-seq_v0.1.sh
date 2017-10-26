#!/bin/bash
set -e
set -u

# Take raw ChIP-seq reads, perform adapter & quality trimming, align to TAIR10 genome
# make sure genome index has been built:
# subread-buildindex -o TAIR10_subread_index TAIR10_chr1.fas TAIR10_chr2.fas TAIR10_chr3.fas TAIR10_chr4.fas TAIR10_chr5.fas TAIR10_chrC.fas TAIR10_chrM.fas

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: chip-seq_v0.1.sh <SE, PE> <fastq R1> <R2> <subread indexed genome> <fileID output>"
exit 1
fi

###
# SINGLE END
###

if [ "$1" == "SE" ]; then

# requirements
if [ "$#" -ne 4 ]; then
echo "Missing required arguments for single-end!"
echo "USAGE: chip-seq_v0.1.sh <SE> <R1> <subread indexed ref genome> <fileID output>"
exit 1
fi

#gather input variables
type=$1
fq=$2;
index=$3; #path to subread indexed reference genome
fileID=$4;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing ChIP-seq alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_ChIP_${dow}
mv $fq ${fileID}_ChIP_${dow}
cd ${fileID}_ChIP_${dow}

if [[ $fq != *.gz ]];then
gzip $fq
fq="${fq}.gz"
fi

# initial fastqc
mkdir 1_fastqc
fastqc $fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq%%.fastq*}_fastqc* 1_fastqc

# adapter and quality trimming using scythe and sickle 
mkdir 2_scythe_sickle
cd 2_scythe_sickle

scythe -a /home/diep/scripts/TruSeq-adapters.fa -p 0.1 ../$fq > ${fq%%.fastq*}_noadapt.fastq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

sickle se -f ${fq%%.fastq*}_noadapt.fastq -o ${fq%%.fastq*}_trimmed.fastq -t sanger -q 20 -l 20 2>&1 | tee -a ../${fileID}_logs_${dow}.log

rm ${fq%%.fastq*}_noadapt.fastq
cd ../

# fastqc again
mkdir 3_trimmed_fastqc
fastqc 2_scythe_sickle/${fq%%.fastq*}_trimmed.fastq 2>&1 | tee -a ${fileID}_logs_${dow}.log

mv 2_scythe_sickle/${fq%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc

mkdir 0_fastq
mv $fq 0_fastq

# subread align
mkdir 4_subread-align
mv 2_scythe_sickle/${fq%%.fastq*}_trimmed.fastq -t 4_subread-align/
cd 4_subread-align/

echo "Beginning alignment ..."

## Subread aligner
subread-align -T 2 -t 1 -u -H -i $index -r ${fq%%.fastq*}_trimmed.fastq -o "${fileID}.bam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

## Subjunc aligner
## subjunc -T 2 -u -H -i $index -r ${fq%%.fastq*}_trimmed.fastq -o "${fileID}.bam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

echo "Cleaning ..."

if [[ $fq%%.fastq}* != *".gz" ]]; then gzip ${fq%%.fastq*}_trimmed.fastq; fi

tmpbam="${fileID}.bam"
outbam="${fileID}.sorted.bam"
samtools sort -m 2G ${tmpbam} -o "${outbam}" 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools index ${outbam}.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
rm -v ${tmpbam}
mv *trimmed.fastq.gz ../2_scythe_sickle/

echo "Alignment complete"

fi

#### 
# PAIRED END
####

if [ "$1" == "PE" ]; then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: chip-seq_v0.1.sh <PE> <R1> <R2> <subread indexed genome> <fileID output>"
exit 1
fi

#gather input variables
type=$1
fq1=$2;
fq2=$3;
index=$4; #path to genome index
fileID=$5;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing ChIP-seq alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq1 $fq2"
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
gzip $fq1
fq1="${fq1}.gz"
fi

if [[ $fq2 != *.gz ]];then
gzip $fq2
fq2="${fq2}.gz"
fi

# initial fastqc
mkdir 1_fastqc
fastqc $fq1 $fq2 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq1%%.fastq*}_fastqc* 1_fastqc
mv ${fq2%%.fastq*}_fastqc* 1_fastqc

# adapter and quality trimming using scythe and sickle 
mkdir 2_scythe_sickle
cd 2_scythe_sickle

scythe -a /home/diep/scripts/TruSeq-adapters.fa -p 0.1 ../$fq1 > ${fq1%%.fastq*}_noadapt.fastq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

scythe -a /home/diep/scripts/TruSeq-adapters.fa -p 0.1 ../$fq2 > ${fq2%%.fastq*}_noadapt.fastq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

sickle pe -f ${fq1%%.fastq*}_noadapt.fastq -r ${fq2%%.fastq*}_noadapt.fastq -o ${fq1%%.fastq*}_trimmed.fastq -p ${fq2%%.fastq*}_trimmed.fastq -s trimmed.singles.fastq -t sanger -q 20 -l 20 2>&1 | tee -a ../${fileID}_logs_${dow}.log

rm ${fq1%%.fastq*}_noadapt.fastq
rm ${fq2%%.fastq*}_noadapt.fastq
cd ../

# fastqc again
mkdir 3_trimmed_fastqc
fastqc 2_scythe_sickle/${fq1%%.fastq*}_trimmed.fastq 2_scythe_sickle/${fq2%%.fastq*}_trimmed.fastq 2>&1 | tee -a ${fileID}_logs_${dow}.log

mv 2_scythe_sickle/${fq1%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc
mv 2_scythe_sickle/${fq2%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc

mkdir 0_fastq
mv $fq1 0_fastq
mv $fq2 0_fastq

# subread align
mkdir 4_subread-align
mv 2_scythe_sickle/${fq1%%.fastq*}_trimmed.fastq -t 4_subread-align/
mv 2_scythe_sickle/${fq2%%.fastq*}_trimmed.fastq -t 4_subread-align/
cd 4_subread-align/

echo "Beginning alignment ..."

tmpbam="${fileID}.bam"
outbam="${fileID}.sorted.bam"

subread-align -T 2 -H -t 1 -u -i ${index} -r ${fq1%%.fastq*}_trimmed.fastq -R ${fq2%%.fastq*}_trimmed.fastq -o "${tmpbam}" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

echo "cleaning..."

if [[ $fq1%%.fastq}* != *".gz" ]]; then gzip ${fq1%%.fastq*}_trimmed.fastq; fi
if [[ $fq2%%.fastq}* != *".gz" ]]; then gzip ${fq2%%.fastq*}_trimmed.fastq; fi

samtools sort -m 2G ${tmpbam} $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools index ${outbam}.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# delete temp bam and gzip sam
rm -v ${tmpbam}
mv *trimmed.fastq.gz ../2_scythe_sickle/

fi
