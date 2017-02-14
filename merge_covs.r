#!/usr/bin/env Rscript

# merge cov files to make pairwise correlation matrices per cytosine base

options(echo=T)
args = commandArgs(trailingOnly=T)
print(args)
context=args[1]

files=dir(pattern=paste0(context,"*.bed.bismark.cov"))

data <- read.delim(files[1], head=F)
data <- data[,1:4]
data <- data[data$V1 != 'Mt',]
data <- data[data$V1 != 'Pt',]
test <- as.numeric(regexec(text=paste0(files[1]), pattern='_C'))
name <- substr(paste0(files[1]), start=1, stop=test-1)
colnames(data)=c('V1','V2','V3',paste(name))

for(i in 2:length(files)){
file=read.delim(files[i],head=F)
file=file[,1:4]
file <- file[file$V1 != 'Mt',]
file <- file[file$V1 != 'Pt',]
test <- as.numeric(regexec(text=paste0(files[i]), pattern='_C'))
name <- substr(paste0(files[i]), start=1, stop=test-1)
colnames(file)=c('V1','V2','V3',paste(name))
temp=merge(data,file,by=c('V1','V2','V3'),all=T)
data=temp
}

test=data[complete.cases(data),]
a <- cor(as.matrix(test[,4:length(test)]))

# write out correlation matrix with samples hc'd
# hc <- hclust(as.dist(a))
# b <- a[hc$order, hc$order]
# write.table(b, 'correlation_matrix_hc_ordered.txt', sep='\t', row.names=T, col.names=T, quote=F)

library(gplots)

pdf(paste0(context,"-bp_cor.pdf"))
heatmap.2(a,
          trace = 'none',
          density.info = "none",
          symm = F,
          symkey = F,
          key = T,
          dendrogram='both',
          )
dev.off()

# PNG output
# png(file=paste0('bp_cor_',context,'.png'), width=800, height = 750, res=300, pointsize = 3)
# heatmap.2(a, trace='none',density.info='none',symm=F,symkey=F,key=T,dendrogram='both',cexCol=1,cexRow=1,srtCol=45)
# dev.off()
