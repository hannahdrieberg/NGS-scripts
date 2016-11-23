#!/bin/bash

# Using TAIR GFF files to produce annotation files
# Derived from SRE gene_to_gene.sh

# Get TAIR annotation file
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff

# Make bedfile file of TAIR10 genes

R

gene_te=read.delim('TAIR10_GFF3_genes_transposons.gff', header=F)
colnames(gene_te) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
gene <- subset(gene_te, feature == 'gene')

gene$chr <- substr(as.character(gene$seqname), start=4, stop=4)
gene$Name <- sapply(strsplit(as.character(gene$attributes), split=';'), function(l) l[3])
test <- as.numeric(regexec(pattern='Name', text=as.character(gene$Name)))
gene$Name <- substr(x=as.character(gene$Name), start = test+5, stop = nchar(as.character(gene$Name)))
gene$ID <- sapply(strsplit(as.character(gene$attributes), split=';'), function(l) l[1])
test <- as.numeric(regexec(pattern='ID', text=as.character(gene$ID)))
gene$ID <- substr(x=as.character(gene$ID), start = test+3, stop = nchar(as.character(gene$ID)))
gene$type <- sapply(strsplit(as.character(gene$attributes), split=';'), function(l) l[2])
test <- as.numeric(regexec(pattern='Note', text=as.character(gene$type)))
gene$type <- substr(x=as.character(gene$type), start = test+5, stop = nchar(as.character(gene$type)))

# Make bedfile file of TAIR10 te

te <- subset(gene_te, feature == 'transposable_element')

te$chr <- substr(as.character(te$seqname), start=4, stop=4)
te$Name <- sapply(strsplit(as.character(te$attributes), split=';'), function(l) l[2])
test <- as.numeric(regexec(pattern='Name', text=as.character(te$Name)))
te$Name <- substr(x=as.character(te$Name), start = test+5, stop = nchar(as.character(te$Name)))
te$ID <- sapply(strsplit(as.character(te$attributes), split=';'), function(l) l[1])
test <- as.numeric(regexec(pattern='ID', text=as.character(te$ID)))
te$ID <- substr(x=as.character(te$ID), start = test+3, stop = nchar(as.character(te$ID)))
te$type <- sapply(strsplit(as.character(te$attributes), split=';'), function(l) l[3])
test <- as.numeric(regexec(pattern='Alias', text=as.character(te$type)))
te$type <- substr(x=as.character(te$type), start = test+6, stop = nchar(as.character(te$type)))

gene.out=gene[,c('chr','start','end','Name','score','strand')]
write.table(gene.out,'TAIR10_genes.bed',sep='\t',row.names=F,col.names=F,quote=F)

te.out=te[,c('chr','start','end','Name','score','strand')]
write.table(te.out,'TAIR10_TE.bed',sep='\t',row.names=F,col.names=F,quote=F)

gene_te.out = rbind(gene.out, te.out)
write.table(gene_te.out,'TAIR10_genes_TE.bed',sep='\t',row.names=F,col.names=F,quote=F)

# make bedfile of TAIR10 5'UTR

futr <- subset(gene_te, feature == 'five_prime_UTR')

futr$chr <- substr(as.character(futr$seqname), start=4, stop=4)
test <- as.numeric(regexec(pattern='Parent', text=as.character(futr$attributes)))
futr$Name <- substr(x=as.character(futr$attributes), start = test+7, stop = nchar(as.character(futr$attributes)))

futr.out=futr[,c('chr','start','end','Name','score','strand')]
write.table(futr.out,'TAIR10_5_utr.bed',sep='\t',row.names=F,col.names=F,quote=F)

quit()
n
