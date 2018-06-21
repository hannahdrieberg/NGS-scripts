#!/usr/bin/env Rscript
options(echo=T)
library(fields)
args=commandArgs(trailingOnly=T)
print(args)

## Citation
# Taudt, A., Roquis, D., Vidalis, A., Wardenaar, R., Johannes, F., and Colome-Tatché́-Tatché, M. (2018). METHimpute: imputation-guided construction of complete methylomes from WGBS data. BMC Genomics 19: 444.

## Perform METHimpute to get imputed/recalibrated genome-wide methylation levels at single Cs and 100bp tiles
# https://github.com/ataudt/methimpute/blob/master/README.md
# https://github.com/ataudt/methimpute/blob/master/vignettes/methimpute.pdf

### Installation
# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("GenomicRanges"))
# library(devtools)
# install_github("ataudt/methimpute")

### Input files
# Run in bash to get 1-based genome-wide cytosine report
# bismark_methylation_extractor --multicore 4 --cytosine_report --CX --genome_folder $HOME/TAIR10 *sorted.bam

## load library
library(methimpute)
library(tidyverse)

## file "CX_report.txt"
file <- args[1]
outname <- args[2]

## chromosome lengths from methimpute
data(arabidopsis_chromosomes)
arabidopsis_chromosomes$chromosome <- sub('chr', 'Chr', arabidopsis_chromosomes$chromosome)

## data import
bismark.data <- importBismark(file, chrom.lengths=arabidopsis_chromosomes)

## Get positions of all cytosines to inflate methylation data (include non-covered sites)
fasta.file <- '~/TAIR10/chromosomes/arabidopsis_seq.fa'
cytosine.positions = extractCytosinesFromFASTA(fasta.file, contexts = c('CG','CHG','CHH'))
methylome = inflateMethylome(bismark.data,cytosine.positions)
print(methylome)

## Obtain correlation parameters (methylation levels from adjacent cytosines)
distcor = distanceCorrelation(methylome, separate.contexts = TRUE)

## Estimate decay parameter for distancce dependeny of the transition probabilities in HMM
fit = estimateTransDist(distcor)

## HMM for complete set using transition probabilities
model = callMethylationSeparate(data = methylome, transDist = fit$transDist, num.threads = 4)
# print(model)

## At genes and TE coordinates
data(arabidopsis_genes)
seqlevels(arabidopsis_genes) <-  sub('chr', 'Chr', seqlevels(arabidopsis_genes))
data(arabidopsis_TEs)
seqlevels(arabidopsis_TEs) <- sub('chr', 'Chr', seqlevels(arabidopsis_TEs))

## METHimpute plotting
pdf(paste0(outname,"_methimpute_HMMfit_enrichment.pdf"))
print(fit$plot)
plotHistogram(model, total.counts=5)
plotScatter(model)
plotTransitionProbs(model)
plotConvergence(model)
plotPosteriorDistance(model$data)
plotEnrichment(model, annotation=arabidopsis_genes)
plotEnrichment(model, annotation=arabidopsis_TEs)
dev.off()

## Export full fitted HMM model
# exportMethylome(model, paste0(outname,"_methimpute_HMMfit.tsv"))

## Output recalibrated methylation levels for downstream analysis akin to bismark cov files
df <- methods::as(model$data, 'data.frame') %>%
select(seqnames, start, end, context, rc.meth.lvl)

df_CG <- subset(df, context == "CG") %>%
select(-context) %>%
mutate(rc.meth.lvl = rc.meth.lvl * 100) %>%
utils::write.table(., file = paste0(outname,"_recal_CpG.bed.cov"), quote = F, sep = '\t', row.names = F, col.names = F)

df_CHG <- subset(df, context == "CHG") %>%
select(-context) %>%
mutate(rc.meth.lvl = rc.meth.lvl * 100) %>%
utils::write.table(., file = paste0(outname,"_recal_CHG.bed.cov"), quote = F, sep = '\t', row.names = F, col.names = F)

df_CHH <- subset(df, context == "CHH") %>%
select(-context) %>%
mutate(rc.meth.lvl = rc.meth.lvl * 100) %>%
utils::write.table(., file = paste0(outname,"_recal_CHH.bed.cov"), quote = F, sep = '\t', row.names = F, col.names = F)

## Binned methylation output of recalibrated weighted methylation levels
df_100bp <- binMethylome(model$data, binsize=100, contexts=c("CG","CHG","CHH"), columns.average="rc.meth.lvl")

df_100bp_CG <- methods::as(df_100bp$CG, 'data.frame') %>%
select(seqnames, start, end, rc.meth.lvl) %>%
mutate(rc.meth.lvl = rc.meth.lvl * 100) %>%
mutate(start = start - 1) %>%
mutate(end = end - 1) %>%
utils::write.table(., file = paste0(outname,"_recal_CpG_100bp.bed"), quote = F, sep = '\t', row.names = F, col.names = F)

df_100bp_CHG <- methods::as(df_100bp$CHG, 'data.frame') %>%
select(seqnames, start, end, rc.meth.lvl) %>%
mutate(rc.meth.lvl = rc.meth.lvl * 100) %>%
mutate(start = start - 1) %>%
mutate(end = end - 1) %>%
utils::write.table(., file = paste0(outname,"_recal_CHG_100bp.bed"), quote = F, sep = '\t', row.names = F, col.names = F)

df_100bp_CHH <- methods::as(df_100bp$CHH, 'data.frame') %>%
select(seqnames, start, end, rc.meth.lvl) %>%
mutate(rc.meth.lvl = rc.meth.lvl * 100) %>%
mutate(start = start - 1) %>%
mutate(end = end - 1) %>%
utils::write.table(., file = paste0(outname,"_recal_CHH_100bp.bed"), quote = F, sep = '\t', row.names = F, col.names = F)

