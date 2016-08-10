#!/usr/bin/env Rscript
# Use wigs to look at methylation across genome in samples 
# Provide context of interest

options(echo=T)
args = commandArgs(trailingOnly=T)
print(args)
context=args[1]

files=dir(pattern=paste0(context,"_100bp.wig"))
data <- read.delim(files[1], head=F, skip=1)
data <- data[,1:4]

colnames(data)=c('V1','V2','V3',paste(files[1]))
for(i in 2:length(files)){
file=read.delim(files[i],head=F,skip=1)
file=file[,1:4]
colnames(file)=c('V1','V2','V3',paste(files[i]))
temp=merge(data,file,by=c('V1','V2','V3'),all=T)
data=temp
}

test=data[complete.cases(data),]
test=test[test$V1 != 'Mt']
test=test[test$V1 != 'Pt']
heat <- test[4:length(test)]
b <- as.numeric(regexec(pattern='_C', text=names(heat))
c <- substr(names(heat), start = 1, stop = b - 1)
names(heat) = c
heat <- as.matrix(heat)
library(gplots)

metpalette = colorRampPalette(c("white","red"))(99)
metbreak = c(seq(0,25,length=25),seq(26,50,length=25),seq(51,75,length=25),seq(76,100,length=25))

rowDistance = dist(heat)
rowCluster = hclust(rowDistance) 
rowDend = as.dendrogram(rowCluster)


heatmap.2(heat,
          col = metpalette,
          breaks = metbreak,
          trace = 'none',
          density.info = "none",
          symm = F,
          symkey = F,
          symbreaks = T,
          scale = "none",
          key = T,
          labRow = F,
          dendrogram='both',
          Colv=T
          )



pdf(paste0("heatmap_"context.pdf"))
heatmap.2(heat,
          col = metpalette,
          breaks = metbreak,
          trace = 'none',
          density.info = "none",
          symm = F,
          symkey = F,
          symbreaks = T,
          scale = "none",
          key = T,
          labRow = F,
          Rowv = rowDend,
          dendrogram='both',
          Colv=T
          )
dev.off()
