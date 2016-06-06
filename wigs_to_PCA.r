#!/usr/bin/env Rscript
# R script to take all wig files in a directory to construct a correlation matrix 
# between samples based on average methylation across 100bp windows. 
# Then perform PCA.
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
a <- cor(as.matrix(test[,4:length(test)]))
dd <- as.dist(a)

pdf(file=paste0(context,"_hclust_tree.pdf"))
plot(hclust(dd))
dev.off()

name = as.numeric(regexec(pattern="_C", text=colnames(a)))
colnames(a) = substr(x=colnames(a), start=1, stop=name-1)

pc=prcomp(a)
pdf(file=paste0(context,"_pca.pdf"))
plot(pc, type ='l' , main='Variance of PCs')
plot(pc$x[1,], pc$x[2,], xlab = 'PC1', ylab='PC2')
text(pc$x[1,], pc$x[2,], colnames(a), cex = 0.8, pos=4)
dev.off()


