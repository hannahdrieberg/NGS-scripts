## Template file to perform DEG calling using for RNAseq data with edgeR

## Installation of edgeR
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('edgeR'))

library(tidyverse)
library(edgeR)

## Count files from featureCounts
countFiles <- dir(pattern = ".counts")
## Define sample groups
sampleGroups <- c(1,1,1,2,2,2)

## The DGElist object
lbls <- sapply(strsplit(countFiles, "_"), function(l) l[1])
dge <- readDGE(countFiles, columns = c(1,7), group = sampleGroups, label=lbls, skip=1)
geneNames <- as.character(rownames(dge$counts))

## rRNA contamination
rRNA <- read.delim("~/scripts/At_rRNA_AGIs.txt", head=F)
rRNA <- as.character(rRNA$V1)
## find rRNAs & count
rRNA.tags <- match(rRNA, geneNames)
rRNA_counts <- dge$counts[rRNA.tags, ]
rRNA.summary <- colSums(rRNA_counts)
rRNA.rates <- (rRNA.summary/dge$samples$lib.size)

# plot(rRNA.rates, main = "rRNA contamination in libraries", ylab = "rRNA abundance (% of total mapped reads)",
xlab = "Sample")
# dev.off()

## rRNA filter
dge$counts <- dge$counts[-rRNA.tags, ]

## Abundance filter (CPM > 1 in at least 3 samples)
keep <- rowSums(cpm(dge) > 1) > 3
dge <- dge[keep, ]

## Re-calculate lib size based on kept transcripts
dge$samples$lib.size <- colSums(dge$counts)

## TMM Normalization
dge.tmm <- calcNormFactors(dge, method = "TMM")

## Estimate common, trended and tagwise dispersion
dge.tmm.disp <- estimateDisp(dge.tmm)

## Use Araport11 annotation file to obtain gene lengths
ara11 <- read.delim("~/Araport11/annotations/Araport11_mRNA.sorted.bed", head=F) %>%
mutate(length = V3 - V2) %>%
filter(V4 %in% geneNames)
gene.lengths <- ara11$length[ara11$V4 %in% rownames(dge.tmm.disp$counts)]

## calculate CPM by group
logcpm <- cpmByGroup(dge.tmm.disp, prior.count=2, log=TRUE, normalized.lib.sizes=TRUE, dispersion=dge.tmm.disp$trended.dispersion)
rpkm_gr <- rpkmByGroup(dge.tmm.disp, gene.length=gene.lengths, dispersion=dge.tmm.disp$trended.dispersion)

## single factor exact tests
et <- exactTest(dge.tmm.disp, pair=levels(dge.tmm.disp$samples$group), dispersion="trended")
# Extract top DEGs based on exactTest
tt <- topTags(et, adjust.method = "fdr", sort.by="logFC", p.value=0.05)
# Identify significantlly differentially expressed genes
sigdeg <- decideTestsDEG(et, adjust.method = "fdr", p.value = 0.05, lfc = 1)
summary(sigdeg)

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

