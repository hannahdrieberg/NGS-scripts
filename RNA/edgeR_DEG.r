## Template file to perform DEG calling using for RNAseq data with edgeR

## Installation of edgeR
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('edgeR'))

library(tidyverse)
library(edgeR)

## Count files from featureCounts
countFiles <- dir(pattern = ".counts")
library(tidyverse)
## Define sample groups
sampleGroups <- c(1,1,1,2,2,2)

## The DGElist object
lbls <- sapply(strsplit(countFiles, "_"), function(l) l[1])
dge <- readDGE(countFiles, columns = c(1,7), group = sampleGroups, label=lbls, skip=1)
geneNames <- as.character(rownames(dge$counts))

## rRNA filter
rRNA <- read.delim("/home/diep/scripts/At_rRNA_AGIs.txt", head=F)
rRNA <- as.character(rRNA$V1)
## find rRNAs & count
rRNA.tags <- match(rRNA, geneNames)
rRNA_counts <- dge$counts[rRNA.tags, ]
rRNA.summary <- colSums(rRNA_counts)
rRNA.rates <- (rRNA.summary/dge$samples$lib.size)

# plot(rRNA.rates, main = "rRNA contamination in libraries", ylab = "rRNA abundance (% of total mapped reads)",
xlab = "Sample")
# dev.off()

dge$counts <- dge$counts[-rRNA.tags, ]

## Abundance filter (CPM > 1 in at least 3 samples)
keep <- rowSums(cpm(dge) > 1) > 3
dge <- dge[keep, ]

## Re-calculate lib size based on kept transcripts
dge$samples$lib.size <- colSums(dge$counts)

## Normalization TMM
dge.tmm <- calcNormFactors(dge, method = "TMM")

## Estimate common, trended and tagwise dispersion
dge.tmm.disp <- estimateDisp(dge.tmm)

## single factor exact tests
et <- exactTest(dge.tmm.disp, pair=1:2)
topTags(et, adjust.method = "fdr", sort.by="logFC", p.value=0.05)
de.tmm <- decideTests(et, adjust.method = "fdr", p.value = 0.05, lfc = 1)
summary(de.tmm)

## Plot biological coefficient of variation and multidimensional scaling plot
# plotBCV(dge.tmm.disp)
# groups <- unique(as.character(sampleGroups))
# n.reps = length(sampleGroups)/length(groups)
# plotMDS(dge.tmm.disp, dim.plot = c(1,2), col = rep(rainbow(length(groups)), each = n.reps))
# plotMDS(dge.tmm.disp, dim.plot = c(2,3), col = rep(rainbow(length(groups)), each = n.reps)) 
# dev.off()

## GLM
## Estimate GLM dispersion paramters

# design <- model.matrix(~0 + sampleGroups)
glm.dge.tmm <- estimateGLMCommonDisp(dge.tmm) 
glm.dge.tmm <- estimateGLMTrendedDisp(glm.dge.tmm)
glm.dge.tmm <- estimateGLMTagwiseDisp(glm.dge.tmm)

# plotBCV(glm.dge.tmm)
# dev.off()

fit <- glmFit(glm.dge.tmm)

