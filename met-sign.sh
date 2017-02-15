#!/bin/bash
set -u

# re-extract DNA methylation, at custom sequence contexts, from sam file and produce cytosine report
# perform in 4_bismark output sub-directory of wgbs workflow

if [ "$#" -ne 2 ]; then
echo "USAGE: met-sign.sh <context> <file> "
exit 1
fi

context=$1
fl=$2

gzip -d ${file}

if [context -eq CHH]; then
seq = (CAA, CAC, CAT, CCA, CCC, CCT, CTA, CTC, CTT)
fi

bismark_methylation_extractor --comprehensive --cytosine_report --CX --genome_folder ~/TAIR10_bs/  --report --buffer_size 10G -s ${fls}

# cleanup
gzip *sam
mkdir seq-contexts
mv *CX_report.txt seq-contexts
rm *txt
cd seq-contexts

# use grep to get output of specific sequence context from report files
for SEQ in seq
do
grep -e $seq -f *.CX_report.txt > out_${seq}.bed
done

