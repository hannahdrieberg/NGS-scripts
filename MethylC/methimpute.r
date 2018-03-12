options(echo=T)
library(fields)
args=commandArgs(trailingOnly=T)
print(args)

## Perform METHimpute
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
# bismark_methylation_extractor --multicore 4 --cytosine_report --CX --genome_folder /home/diep/TAIR10 *sorted.bam

## load library
library(methimpute)

## file "CX_report.txt"
file <- args[1]

## chromosome lengths from methimpute
data(arabidopsis_chromosomes)
arabidopsis_chromosomes$chromosome <- sub('chr', 'Chr', arabidopsis_chromosomes$chromosome)

## data import
bismark.data <- importBismark(file, chrom.lengths=arabidopsis_chromosomes)

## Get positions of all cytosines to inflate methylation data (include non-covered sites)
fasta.file <- '/home/diep/TAIR10/arabidopsis_seq.fa'
cytosine.positions = extractCytosinesFromFASTA(fasta.file, contexts = c('CG','CHG','CHH'))
methylome = inflateMethylome(bismark.data,cytosine.positions)
print(methylome)

## obtain correlation parameters (methylation levels from adjacent cytosines)
distcor = distanceCorrelation(methylome)
# estimate decay parameter for distancce dependeny of the transition probabilities in HMM
fit = estimateTransDist(distcor)
# print(fit)
# dev.off()

## HMM for complete set using transition probabilities
model = callMethylation(data = methylome, transDist = fit$transDist)
print(model)

## Per context HMMs
# model <- callMethylationSeparate(data=methylome)
# print(model)

### METHimpute plotting
outname <- sapply(strsplit(file, "_"), function(l) l[1])
pdf(paste0(outname"_methimpute_HMMfit.pdf"))
print(fit)
plotHistogram(model, total.counts=5)
plotScatter(model)
plotTransitionProbs(model)
plotConvergence(model)
plotPosteriorDistance(model$data)
dev.off()

# At genes and TE coordinates
data(arabidopsis_genes)
seqlevels(arabidopsis_genes) <-  sub('chr', 'Chr', seqlevels(arabidopsis_genes))
data(arabidopsis_TEs)
seqlevels(arabidopsis_TEs) <- sub('chr', 'Chr', seqlevels(arabidopsis_TEs))

# At enrichment plots
pdf("METHimpute_enrichments_plots.pdf")
plotEnrichment(model, annotation=arabidopsis_genes)
plotEnrichment(model, annotation=arabidopsis_TEs)
dev.off()

## Bin cytosine positions, methylation counts in equidistant windows
# binMethylome(data, binsize, contexts = "total", columns.average = c("rc.meth.lvl"))

## Output recalibrated methylation levels
exportMethylome(model, filename = paste0(outname,".tsv"))
CG_out <- model$data[model$data$context == "CG"]
CHG_out <- model$data[model$data$context == "CHG"]
CHH_out <- model$data[model$data$context == "CHH"]

