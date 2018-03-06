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
# bismark_methylation_extractor --cytosine_report --CX --genome_folder /home/diep/TAIR10 *sorted.bam

## load library
library(methimpute)

## actual files
file <- dir(pattern="*.CX_report.txt")

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
fit = estimateTransDist(distcor)

print(fit)
dev.off()

model = callMethylation(data = methylome, transDist = fit$transDist)
plot(model)

