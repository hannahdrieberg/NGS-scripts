#!/bin/bash
set -u
# execute R file DSS_file_prep.sh to make files to run through DMR caller: DSS
# run in folder with cov files of interest and tell which context you want to merge together


if [ "$#" -lt 1 ]; then
echo "### USE = DSS_file_prep.sh {context} ###"
exit 1
fi

context=$1

Rscript /home/diep/scripts/DSS_file_prep.r ${context}

mkdir DSS_${context}_output
mv *_${context}_output.txt DSS_${context}_output/
