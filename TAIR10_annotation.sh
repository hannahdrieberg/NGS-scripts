#!/bin/bash

# Use ENSEMBL TAIR GFF files to produce annotation files
wget ftp://ftp.ensemblgenomes.org/pub/release-39/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.39.gff3.gz

# Make bedfile file of TAIR10 genes

R

library(tidyverse)

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

ens <- gffRead('Arabidopsis_thaliana.TAIR10.39.gff3') %>%
	mutate(seqname=sub("Mt","M", seqname)) %>%
	mutate(seqname=sub("Pt","C", seqname)) %>%
	mutate(seqname=paste0('Chr',seqname))

# Chromosome annotation
chr <- subset(ens,ens$feature=='chromosome') %>%
        select(c('seqname','start','end'))

# Gene annotation
gene <- subset(ens, ens$feature == 'gene') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(gene,'TAIR10.39_genes.bed', sep='\t', row.names=F, col.names=F, quote=F)

# mRNA annotation
mRNA <- subset(ens, ens$feature == 'mRNA') %>%
	mutate(ID=getAttributeField(attributes, 'ID')) %>%
	mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
	select(c('seqname','start','end','Name','score','strand'))

write.table(mRNA,'TAIR10.39_mRNA.bed', sep='\t', row.names=F, col.names=F, quote=F)

