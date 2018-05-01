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

# Make descriptive group labels
sampleGroups <- sub(1, "unstressed", x=sampleGroups)
sampleGroups <- sub(2, "drought", x=sampleGroups)

## The DGElist object
lbls <- sapply(strsplit(countFiles, "_"), function(l) l[1])
dge <- readDGE(countFiles, columns = c(1,7), group = sampleGroups, label=lbls, skip=1)
geneNames <- as.character(rownames(dge$counts))

# Use sample groups to make design matrix
design <- model.matrix(~0 + sampleGroups)
colnames(design) <- unique(sampleGroups)

## rRNA contamination
rRNA <- read.delim("~/scripts/At_rRNA_AGIs.txt", head=F)
rRNA <- as.character(rRNA$V1)
## find rRNAs & count
rRNA.tags <- match(rRNA, geneNames)
rRNA_counts <- dge$counts[rRNA.tags, ]
rRNA.summary <- colSums(rRNA_counts)
rRNA.rates <- (rRNA.summary/dge$samples$lib.size)

## rRNA filter
dge$counts <- dge$counts[-rRNA.tags, ]

## Abundance filter (CPM > 1 in at least 3 samples)
keep <- rowSums(cpm(dge) > 1) > 3
dge <- dge[keep, ]

## Re-calculate lib size based on retained transcripts
dge$samples$lib.size <- colSums(dge$counts)

## TMM Normalization
dge.tmm <- calcNormFactors(dge, method = "TMM")

##########
## exact tests
##########

## Estimate common, trended and tagwise dispersion
dge.tmm.disp <- estimateDisp(dge.tmm, design)

# diagnostic plots
pdf("exacttests_diagnostic_plots.pdf", paper="a4r")
plot(rRNA.rates, main = "rRNA contamination in libraries", ylab = "rRNA abundance (% of total mapped reads)",
xlab = "Sample")
plotBCV(dge.tmm.disp)
groups <- unique(as.character(sampleGroups))
n.reps <- length(sampleGroups)/length(groups)
plotMDS(dge.tmm.disp, dim.plot = c(1,2), col = rep(rainbow(length(groups)), each = n.reps))
plotMDS(dge.tmm.disp, dim.plot = c(2,3), col = rep(rainbow(length(groups)), each = n.reps))
dev.off()

## Use Araport11 annotation file to obtain gene lengths
ara11 <- read.delim("~/Araport11/annotations/Araport11_mRNA.sorted.bed", head=F) %>%
mutate(length = V3 - V2) %>%
filter(V4 %in% geneNames)
gene.lengths <- ara11$length[ara11$V4 %in% rownames(dge.tmm.disp$counts)]

## calculate CPM by group
logcpm <- cpmByGroup(dge.tmm.disp, prior.count=2, log=TRUE, normalized.lib.sizes=TRUE, dispersion=dge.tmm.disp$trended.dispersion)
rpkm_gr <- rpkmByGroup(dge.tmm.disp, gene.length=gene.lengths, dispersion=dge.tmm.disp$trended.dispersion)

## single factor exact tests (pairwise comparisons)
et <- exactTest(dge.tmm.disp, pair=as.character(unique(dge.tmm.disp$samples$group)), dispersion="trended")

# Identify significantlly differentially expressed genes
sigdeg <- decideTestsDGE(et, adjust.method = "fdr", p.value = 0.05, lfc = 1)
summary(sigdeg)

## final DEG (from pairwise et)
# Extract top DEGs based on exactTest
# tt <- topTags(et, adjust.method = "fdr", sort.by="logFC", p.value=0.05)
tt <- topTags(et, adjust.method = "fdr", sort.by="logFC", p.value=0.05, n=dim(et)[1])
tt <- tt$table[abs(tt$table$logFC) >= 1,]

## Output table

getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
        cat("found", nrow(gff), "rows with classes:",
        paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

# read in Araport11 gff3
ara=gffRead("~/Araport11/Araport11_GFF3_genes_transposons.201606.gff")

# Gene annotation
gene=subset(ara,ara$feature=='gene') %>%
	mutate(locus=getAttributeField(attributes, 'ID')) %>%
	mutate(description=getAttributeField(attributes, 'Note')) %>%
	mutate(type=getAttributeField(attributes, 'locus_type')) %>%
	mutate(primary=getAttributeField(attributes, 'full_name')) %>%
	mutate(alias=getAttributeField(attributes, 'Alias')) %>%
	select('locus','description','type','primary','alias') %>%
	subset(locus %in% rownames(tt)) %>%
	mutate(logFC=tt$logFC[match(locus, rownames(tt))]) %>%
	mutate(logCPM=tt$logCPM[match(locus, rownames(tt))]) %>%
	mutate(PVal=tt$PValue[match(locus, rownames(tt))]) %>%
	mutate(FDR=tt$FDR[match(locus, rownames(tt))]) %>%
	arrange(desc(logFC))

write.table(gene,paste0('RNAseq_edgeR-exacttest_',groups[1],'vs',groups[2],'.txt'), sep='\t', na='', col.names=T, row.names=F, quote=F)

################
## GLM
###############

### Estimate GLM dispersion paramters
glm.dge.tmm <- estimateGLMCommonDisp(dge.tmm)
glm.dge.tmm <- estimateGLMTrendedDisp(glm.dge.tmm)
glm.dge.tmm <- estimateGLMTagwiseDisp(glm.dge.tmm)

# diagnostic plots
pdf("exacttests_diagnostic_plots.pdf", paper="a4r")
plot(rRNA.rates, main = "rRNA contamination in libraries", ylab = "rRNA abundance (% of total mapped reads)",
xlab = "Sample")
plotBCV(glm.dge.tmm)
groups <- unique(as.character(sampleGroups))
n.reps <- length(sampleGroups)/length(groups)
plotMDS(dge.tmm.disp, dim.plot = c(1,2), col = rep(rainbow(length(groups)), each = n.reps))
plotMDS(dge.tmm.disp, dim.plot = c(2,3), col = rep(rainbow(length(groups)), each = n.reps))
dev.off()

fit <- glmFit(glm.dge.tmm)


