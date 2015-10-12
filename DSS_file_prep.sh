#!/bin/bash

# execute R file DSS_file_prep.sh to make files to run through DMR caller: DSS
# run in folder with cov files of interest and tell which context you want to merge together

context=$1

Rscript /home/diep/scripts/DSS_file_prep.r ${context}


