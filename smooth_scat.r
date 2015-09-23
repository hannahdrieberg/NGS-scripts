# R script for scatman_smooth.sh to produce smoothed scatterplots for methylation levels vs feature length
# Will produce plots in the order: CpG, CHG, CHH

args=commandArgs(trailingOnly=T)
print(args)

#cpg
a <- read.delim(paste0(args[1],"_CpG_",args[2],".bed"), header=F)
a$length <- a$V7 - a$V6

#chg
b <- read.delim(paste0(args[1],"_CHG_",args[2],".bed"), header=F)
b$length <- b$V7 - b$V6

#chh
c <- read.delim(paste0(args[1],"_CHH_",args[2],".bed"), header=F)
c$length <- c$V7 - c$V6

pdf(file=paste0(args[1],"_",args[2],".pdf"))
par(mfrow=c(2,2))
smoothScatter(x=a$length, y=a$V4, ylab="CpG Methylation", xlab="TE Length (bp)")
smoothScatter(x=b$length, y=b$V4, ylab="CHG Methylation", xlab="TE Length (bp)")
smoothScatter(x=c$length, y=c$V4, ylab="CHH Methylation", xlab="TE Length (bp)")
dev.off()

