#!/usr/bin/env Rscript
# merge wigs to make correlation matrices
# run in directory with wig files to correlate

args = commandArgs(trailingOnly=T)
print(args)
context=args[1]

library(tidyverse)
library(gplots)

files <- dir(pattern=paste0(context,"_100bp.wig"))
data <- data_frame(files) %>%
mutate(file_contents = map(files, read_delim, delim='\t', col_names=F, skip=1)) %>%
unnest() %>%
filter(X1 != 'Mt' & X1 != 'ChrM' & X1 != 'Pt' & X1 != 'ChrC') %>%
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
cor() %>%
as.matrix()

pdf(file=paste0('wig_cor_',context,'.pdf'), width = 0, height = 0, paper="a4")
heatmap.2(data, 
	trace='none',
	density.info='none',
	symm=F,
	symkey=F,
	key=T,
	colsep = 1:ncol(data),
	rowsep = 1:nrow(data),
	sepcolor = "white",
	sepwidth = c(0.001,0.001),
	dendrogram='both')
dev.off()
