# produce mean 5mC levels for R plotting
options(echo=T)
library(fields)
args=commandArgs(trailingOnly=T)
print(args)

smplname <- as.character(paste0(args[1]))
outname <- as.character(paste0(args[2]))
context <- as.character(paste0(args[3]))

data <- dir(pattern=paste0(smplname,"-",outname,"-sub",context,"-report.1k.bed")) %>%
read_delim(delim = '\t', col_names=F) %>%
mutate(rel.dist=ifelse(X13==0,ifelse(X12=="-",((X9-X2)/(X9-X8))*1000,((X2-X8)/(X9-X8))*1000),ifelse(X13>0,X13+1000,X13))) %>%
mutate(fixy=ifelse(rel.dist<0 & X13==0,0,ifelse(rel.dist>1000 & X13==0, 1000, rel.dist)))

out <- NULL
for(i in unique(data$X5)){
a <- subset(data, X5 == i)
a <- stats.bin(a$fixy,a$X6,N=100)
temp <- as.data.frame(cbind(matrix(a$centers,ncol=1),a$stats["mean",]))
temp$motiff <- paste0(i)
out <- rbind(temp, out)
}

write.table(out,paste0(paste(smplname,context,outname,sep='_'),'.txt'),quote=F, col.names=T, row.names=F, sep='\t')
