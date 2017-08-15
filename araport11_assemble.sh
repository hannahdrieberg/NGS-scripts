# build annotations based off gff3 files from araport11 release
# https://www.araport.org/data/araport11
# derived from SRE gene_to_gene.sh
# Run lines manually in annotation directory

# Readme file
wget https://www.araport.org/download_file/Araport11_Release_201606/annotation/README.201606.md

# Araport11 annotation in GFF3

# curl -sO -H 'Authorization: Bearer 745dd29759980b058db8fb9efc7af5' https://api.araport.org/files/v2/media/system/araport-public-files//Araport11_Release_201606/annotation/Araport11_GFF3_genes_transposons.201606.gff.gz

wget http://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz

gzip -d *gff.gz

# Make bed files

R

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

ara=gffRead('Araport11_GFF3_genes_transposons.201606.gff')

# Gene annotation
gene=subset(ara,ara$feature=='gene')
gene$Name=getAttributeField(gene$attributes, 'Name')
gene$ID=getAttributeField(gene$attributes, 'ID')
gene$chr <- sapply(strsplit(gene$seqname, 'Chr'), function(l) l[2])
gene.out=gene[,c('chr','start','end','Name','score','strand')]

write.table(gene.out,'Araport11_genes.bed',sep='\t',row.names=F,col.names=F,quote=F)

# TE annotation
te=subset(ara,ara$feature=='transposable_element')
te$Name=getAttributeField(te$attributes, 'Name')
te$ID=getAttributeField(te$attributes, 'ID')
te$chr <- sapply(strsplit(te$seqname, 'Chr'), function(l) l[2])
te.out=te[,c('chr','start','end','Name','score','strand')]

write.table(te.out,'Araport11_TE.bed',sep='\t',row.names=F,col.names=F,quote=F)

quit()
n

##########

rm *gff

