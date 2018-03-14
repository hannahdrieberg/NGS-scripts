## Template file to perform DEG calling using for RNAseq data with edgeR

## Installation of edgeR
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('edgeR'))

library(edgeR)
library(tidyverse)

## Count files from featureCounts
countFiles <- dir(pattern = ".counts")
## Define sample groups
sampleGroups <- c(1,1,1,2,2,2)

dge <- readDGE(countfiles,

# dge <- readDGE(countFiles, columns = c(1, 3), group = sampleGroups,
# labels = as.character(samples))
# geneNames <- as.character(rownames(dge$counts))


