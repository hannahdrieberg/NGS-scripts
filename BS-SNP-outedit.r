fl = dir()
for(i in 1:length(fl)){
a <- read.delim(fl[i], skip=1)
colnames(a)[1] = "#CHROM"
a <- a[complete.cases(a),]
a <- a[1:7]
write.table(a, file=paste0(fl[i],".edit"), sep='\t', quote=F, row.names=F)
}
.
