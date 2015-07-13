#!/bin/bash

# bash script to take bed files, subsetted at features of interest, or across entire chromosomes or genome-wide, and plot heatmaps. This was getting hard to do on the laptop ...
# all this script is doing is giving the files of interest to be taken into R with the appropriate options/sample id's.


if [ "$#" -ne 1 ]; then
echo "#######################################################################################################"
echo "Run in directory containing bed files of interest --> output from methylation_tiling.sh"
echo "BED filenames must have structure: _<context>_DMRs_<annotation>_<subset>.bed"
echo "USAGE: heatmap.sh <annotation name>"
echo "Note: annotation name should be the output name given to bed files from methylation_tiling.sh"
echo "Example: heatmap.sh TE"
echo "Change R script depending on experiment/samples --> treatment defined here: line 37"
echo "Check R script also for output heatmap file names"
echo "#######################################################################################################"
exit 1
fi

annotation=$1

Rscript $HOME/scripts/Exp277_averaging_replicates.r ${annotation}


