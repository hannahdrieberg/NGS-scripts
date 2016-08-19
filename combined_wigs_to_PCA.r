#!/usr/bin/env Rscript
# R script similar to wigs_to_PCA, except taking all contexts at once.
# Then perform PCA.

options(echo=T)

files=dir(pattern=paste0("_100bp.wig"))
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

# Correlation matrix
a <- cor(as.matrix(test[,4:length(test)]))

name = as.numeric(regexec(pattern="_100", text=colnames(a)))
colnames(a) = substr(x=colnames(a), start=1, stop=name-1)

# Use correlation as distance
dd <- as.dist(1-a)
names(dd) = colnames(a)

# PCA 
pc=prcomp(a)

pdf(file=paste0("100bp_combined_hclust_PCA.pdf"))
plot(hclust(dd), main='hclust on 100bp wigs', cex=0.5)
plot(pc, type ='l' , main='Variance of PCs')
plot(pc$x[1,], pc$x[2,], xlab = 'PC1', ylab='PC2', main = 'PCA on 100bp wigs', cex=0.5)
text(pc$x[1,], pc$x[2,], colnames(a), cex = 0.5, pos=4)
dev.off()
