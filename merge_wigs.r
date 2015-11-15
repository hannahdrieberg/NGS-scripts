# merge wigs to make pairwise correlation matrices

files=dir(pattern="*.wig")
data <- read.delim(files[1], head=F, skip=1)
data <- data[,1:4]
colnames(data)=c('V1','V2','V3',paste(files[1],'_prop',sep=''))
for(i in 2:length(files)){
file=read.delim(files[i],head=F,skip=1)
file=file[,1:4]
colnames(file)=c('V1','V2','V3',paste(files[i],'_prop',sep=''))
temp=merge(data,file,by=c('V1','V2','V3'),all=T)
data=temp
}
# make all prop_met columns numeric
for(i in 4:length(data)){
data[,i] <- as.numeric(data[,i])
}

# psuedo correlation matrix

data=read.delim('output.txt',head=T)
cg=data[grep("*CpG",data$contrast),]
a=strsplit(as.character(cg$contrast),"vs")
a=unlist(a)
a=matrix(a,ncol=2,byrow=T)

b=table(a[,1],a[,2])
write.table(b,'test.txt',sep='\t',row.names=T,col.names=T)

test=data[complete.cases(data),]
write.table(cor(as.matrix(test[,4:length(test)])), 'test.txt', sep='\t', row.names=T, col.names=T, quote=F)
