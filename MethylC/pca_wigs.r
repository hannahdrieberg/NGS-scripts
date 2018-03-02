#!/usr/bin/env Rscript
# merge wigs and perform PCAs

args = commandArgs(trailingOnly=T)
print(args)
context=args[1]
method <- paste0(ifelse(args[2] != "pearson" & args[2] != "kendall" & args[2] != "spearman", yes="pearson", no=args[2]))
method <- ifelse(method == "NA", yes="pearson", no=args[2])
print(paste("cor method = " ,method))

library(tidyverse)

files <- dir(pattern=paste0(context,"_100bp.wig"))
data <- data_frame(files) %>%
mutate(file_contents = map(files, read_delim, delim='\t', col_names=F, skip=1)) %>%
unnest() %>%
filter(X1!='Mt'&X1!='ChrM'&X1!='Pt'&X1!='ChrC') %>%
mutate(sample=sapply(strsplit(files, '_'), function(l) l[1])) %>%
mutate(genotype=sapply(strsplit(sample, '-'), function(l) l[1])) %>%
mutate(rep=sapply(strsplit(sample, '-'), function(l) l[2])) %>%
mutate(X1=ifelse(substr(X1, start=1, stop=3)=="Chr",paste0(X1),paste0("Chr",X1))) %>%
na.omit() %>%
group_by(X1, X2, X3, genotype) %>%
summarise(met = mean(X4)) %>%
spread(key=genotype, value=met) %>%
na.omit() %>%
ungroup() %>%
select(-X1, -X2, -X3) %>%
as.matrix()

# PCA analysis
# pc=prcomp(data)
# plot(pc, type ='l' , main='Variance of PCs')
# plot(pc$x[1,], pc$x[2,], xlab = 'PC1', ylab='PC2')
# text(pc$x[1,], pc$x[2,], colnames(data), cex = 0.8, pos=4)
# library(devtools)
# install_github("ggbiplot","vqv")
# library(ggbiplot)
# ggbiplot(pc, obs.scale=1, var.scale=1, groups=ir.species, ellipse = TRUE, circle = TRUE) + theme(legend.direction = 'horizontal', legend.position = 'top')



