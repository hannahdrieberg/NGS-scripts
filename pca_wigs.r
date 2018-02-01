#!/usr/bin/env Rscript
# merge wigs and perform PCAs

args = commandArgs(trailingOnly=T)
context=args[1]
method <- ifelse(args[2] == '', yes="pearson", no=args[2])
print(args)

# PCA analysis
# pc=prcomp(data)
# plot(pc, type ='l' , main='Variance of PCs')
# plot(pc$x[1,], pc$x[2,], xlab = 'PC1', ylab='PC2')
# text(pc$x[1,], pc$x[2,], colnames(data), cex = 0.8, pos=4)
# library(devtools)
# install_github("ggbiplot","vqv")
# library(ggbiplot)
# ggbiplot(pc, obs.scale=1, var.scale=1, groups=ir.species, ellipse = TRUE, circle = TRUE) + theme(legend.direction = 'horizontal', legend.position = 'top')



