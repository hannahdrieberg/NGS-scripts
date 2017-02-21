# produce mean 5mC levels for R plotting
options(echo=T)
library(fields)
args=commandArgs(trailingOnly=T)
print(args)

smplname <- args[1]
outname <- args[2]
context <- args[3]

if(context=="CHH"){seq=list("CAA", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT")}
out <- NULL
for(i in seq){
# read in files
a <- read.delim(paste0(smplname,'-',outname,'-',context,'-',i,'.1k.bed'), head=F)
a <- subset(a,a$V1!='Mt' & a$V1!='Pt')
a <- a[c(1:3, 7,8,12,14:20)]
real.dist=matrix(ifelse(a$V19=='+',-1*a$V20,a$V20),ncol=1)
a=cbind(a,real.dist)
rel.dist=matrix(ifelse(a$real.dist==0,ifelse(a$V19=="-",((a$V16 - a$V2)/(a$V16 - a$V15))*1000,((a$V2 - a$V15)/(a$V16 - a$V15))*1000),ifelse(a$real.dist>0,a$real.dist + 1000,a$real.dist)),ncol=1)
a=cbind(a,rel.dist)
fixy=ifelse(a$rel.dist < 0 & a$real.dist==0,0,ifelse(a$rel.dist >1000 & a$real.dist==0,1000,a$rel.dist))
a$rel.dist=fixy
a.bin=stats.bin(a$rel.dist,a$V12,N=100)
p.a.bin=cbind(matrix(a.bin$centers,ncol=1),a.bin$stats["mean",])
test <- as.data.frame(p.a.bin)
colnames(test)[1] <- 'pos'
test$sample <- paste0(smplname)
test$context <- paste0(context)
test$motiff <- paste0(i)
out <- rbind(out,test)
}
write.table(out,paste(smplname,'_',context,'_',outname,'.txt',sep=''),quote=F, col.names=T, row.names=F, sep='\t')
