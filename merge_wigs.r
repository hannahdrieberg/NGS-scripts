#!/usr/bin/env Rscript
# merge wigs to make correlation matrices
# run in directory with wig files to correlate

options(echo=T)
args = commandArgs(trailingOnly=T)
print(args)
context=args[1]

require(reshape2)
require(gplots)

files=dir(pattern=paste0(context,"_100bp.wig"))

out <- NULL

for(i in 1:length(files)){
data <- read.delim(files[i],head=F,skip=1)
data <- data[,1:4]
data <- data[data$V1 != 'Mt' & data != "ChrM" & data != "M",]
data <- data[data$V1 != 'Pt' & data != "ChrC" & data != "C",]
data <- data[complete.cases(data),]
data$sample <- sapply(strsplit(files[i], '_'), function(l) l[1])
data$context <- sapply(strsplit(files[i], '_'), function(l) l[2])
data$V1 <- ifelse(substr(data$V1, start=1, stop=3) == "Chr", paste0(data$V1),paste0("Chr",data$V1))
out <- rbind(out,data)
}

test <- dcast(out, formula=V1+V2+V3~sample, value.var = "V4")
test <- test[complete.cases(test),]
a <- cor(as.matrix(test[,4:length(test)]))

pdf(file=paste0('wig_cor_',context,'.pdf'), width=8, height = 7.5, pointsize = 10)
heatmap.2(a, 
	trace='none',
	density.info='none',
	symm=F,
	symkey=F,
	key=T,
	colsep = 1:ncol(a),
	rowsep = 1:nrow(a),
	sepcolor = "white",
	sepwidth = c(0.001,0.001),
	dendrogram='both',
	cexCol=0.8,
	cexRow=0.8)
dev.off()
