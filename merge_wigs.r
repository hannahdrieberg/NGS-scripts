#!/usr/bin/env Rscript

# merge wigs to make pairwise correlation matrices
# run in directory with wig files to correlate

files=dir(pattern="*.wig")
data <- read.delim(files[1], head=F, skip=1)
data <- data[,1:4]
data <- data[data$V1 != 'Mt',]
data <- data[data$V1 != 'Pt',]
test <- as.numeric(regexec(text=paste0(files[1]), pattern='_C'))
name <- substr(paste0(files[1]), start=1, stop=test-1)
colnames(data)=c('V1','V2','V3',paste(name))
for(i in 2:length(files)){
file=read.delim(files[i],head=F,skip=1)
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
hc <- hclust(as.dist(a))
write.table(a, 'correlation_matrix.txt', sep='\t', row.names=T, col.names=T, quote=F)

b <- a[hc$order, hc$order]
write.table(b, 'correlation_matrix_hc_ordered.txt', sep='\t', row.names=T, col.names=T, quote=F)

pdf("cor.pdf")
library(gplots)
heatmap.2(a,
          trace = 'none',
          density.info = "none",
          symm = F,
          symkey = F,
          key = T,
          dendrogram='both',
          )
dev.off()

