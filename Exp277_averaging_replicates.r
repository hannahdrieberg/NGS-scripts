# R script for Exp277 Control vs Drought Samples
# Mapped interesting positions (e.g. DMRs) to genomic features = loci that are near/contain DMRs (any context)
# Extracted methylation information within these loci, and upstream, and downstream; using wig files.
# Read this data in and average replicates together to create a data frame to start visualising
# Extracting sample information based on file naming - THIS IS SUBJECT TO CHANGE EVERY EXPERIMENT
# Set working directory to where files are 
# make sure to have files of a particular feature there or make sure naming is appropriate

# R script for Exp277 Control vs Drought Samples
# Mapped interesting positions (e.g. DMRs) to genomic features = loci that are near/contain DMRs (any context)
# Extracted methylation information within these loci, and upstream, and downstream; using wig files.
# Read this data in and average replicates together to create a data frame to start visualising
# Extracting sample information based on file naming - THIS IS SUBJECT TO CHANGE EVERY EXPERIMENT
# Set working directory to where files are 
# make sure to have files of a particular feature there or make sure naming is appropriate

args=commandArgs(trailingOnly=T)
print(args)

library(plyr)
library(reshape2)
library(limma)
library(gplots)
library(fields)
library(RColorBrewer)

p <- substr(args, start = 0, stop = nchar(args[1]))

a <- dir(pattern = p)
out <- NULL
for(i in 1:length(a)){
input=read.delim(a[i], head=F)
input=subset(input,input$V1!='Mt' & 
                     input$V1!='chrMt' & 
                     input$V1!='Pt' & 
                     input$V1!='chrPt') 
real.dist=matrix(ifelse(input[,14]=='+',-1*input[,15],input[,15]),ncol=1)
real.dist=round(real.dist,digits=-2)
input=cbind(input,real.dist)
y <- as.numeric(regexec(pattern = "_C", text = a[i]))
z <- as.numeric(regexec(pattern = "_D", text = a[i]))
input$sample <- substr(x = as.character(a[i]), start = 0, stop = y-1) 
input$context <- substr(x = as.character(a[i]), start = y+1, stop = y+3)
out <- rbind(input, out)
}
########################################################################
out$Treatment <- factor(ifelse(test = out$sample == "D7" | 
                          out$sample == "D8" |
                          out$sample == "D9", "Drought", "Control"))
########################################################################
# average reps
out.reps <- ddply(out, .(V1,V2,V3,V9,V10,V11,V12,V13,V14,V15,real.dist,Treatment,context), summarise,
             Prop_met = mean(V4),
             Coverage = mean(V7),
             C_count = mean(V8),
             n = sum(!is.na(V4)) # a few funny things happening
             )
# Calculate methylation across Annotation (e.g. Gene Body)
out.body <- subset(out.reps,out.reps$real.dist == 0)
rel.dist=matrix(ifelse(out.body[,9]=="-",((out.body[,6]-(out.body[,2]))/(out.body[,6]-out.body[,5]))*1000,(((out.body[,2])-out.body[,5])/(out.body[,6]-out.body[,5]))*1000),ncol=1)
out.body=cbind(out.body,rel.dist)
fixy=ifelse(out.body$rel.dist < 0,0,ifelse(out.body$rel.dist>1000,1000,out.body$rel.dist))
fixy <- round(fixy, digits = -2)
out.body$rel.dist=fixy
out.v1<-out.body
out.v2=dcast(out.body, V12 ~ factor(rel.dist)+Treatment+context, mean, value.var = "Prop_met")

# Heatmaps of CpG, CHH and CHG methylation within loci body
for(w in unique(out.body$Treatment)){
  e <- out.body[out.body$Treatment == w,]
  pdf(file = paste("./Heatmap_",paste(args[1]),"_",w,".pdf"))
  for(i in unique(e$context)){
    q <- dcast(e[e$context == i,], V12 ~ factor(rel.dist), mean, value.var = "Prop_met")
    rownames(q) <- q[,1]
    q <- q[2:length(q)]
    w <- as.matrix(q)
    w[is.na(w)] <- -10
    heatmap.2(w,
              dendrogram="none", 
              Rowv=NULL,
              Colv=NULL,
              trace='none',
              key=T,
              main = paste(i),
              labRow=F, 
              breaks=c(seq(-1,0,length=1),
                       seq(1,10,length=10),
                       seq(11,20,length=10),
                       seq(21,30,length=10),
                       seq(31,40,length=10),
                       seq(41,50,length=10),
                       seq(51,60,length=10),
                       seq(61,70,length=10),
                       seq(71,80,length=10),
                       seq(81,90,length=10),
                       seq(91,100,length=10)),
              symm=F,
              symkey=F,
              symbreaks = T,
              cexCol = 1,
              col=colorRampPalette(c("white","blue"))(100))
  }
  dev.off()
}

# Methylation upstream and downstream of annotations

# Upstream of annotated loci
out.upstream <- subset(out.reps,out.reps$real.dist < 0)

for(w in unique(out.upstream$Treatment)){
e <- out.upstream[out.upstream$Treatment == w,]
pdf(file = paste("./Heatmap_",paste(args[1]),"_",w,".pdf"))
  for(i in unique(e$context)){
  q <- dcast(e[e$context == i,], V12 ~ factor(real.dist), mean, value.var = "Prop_met")
  rownames(q) <- q[,1]
  q <- q[2:length(q)]
  w <- as.matrix(q)
  w[is.na(w)] <- -10
    heatmap.2(w,
          dendrogram="none", 
          Rowv=NULL,
          Colv=NULL,
          trace='none',
          key=T,
          main = paste0(i),
          labRow=F, 
          breaks=c(seq(-1,0,length=1),
                    seq(1,10,length=10),
                   seq(11,20,length=10),
                   seq(21,30,length=10),
                   seq(31,40,length=10),
                   seq(41,50,length=10),
                   seq(51,60,length=10),
                   seq(61,70,length=10),
                   seq(71,80,length=10),
                   seq(81,90,length=10),
                   seq(91,100,length=10)),
          symm=F,
          symkey=F,
          symbreaks = T,
          cexCol = 1,
          col=colorRampPalette(c("white",
                                  "red"))(100))
  }
  dev.off()
}

# Downstream of genes mapped to DMRs
out.downstream <- subset(out.reps,out.reps$real.dist > 0)

for(w in unique(out.downstream$Treatment)){
e <- out.downstream[out.downstream$Treatment == w,]
pdf(file = paste("Figures/Heatmap_",paste(args[1]),"_",w,".pdf"))
  for(i in unique(e$context)){
  q <- dcast(e[e$context == i,], V12 ~ factor(real.dist), mean, value.var = "Prop_met")
  rownames(q) <- q[,1]
  q <- q[2:length(q)]
  w <- as.matrix(q)
  w[is.na(w)] <- -10
    heatmap.2(w,
          dendrogram="none", 
          Rowv=NULL,
          Colv=NULL,
          trace='none',
          key=T,
          main = paste0(i),
          labRow=F, 
          breaks=c(seq(-1,0,length=1),
                    seq(1,10,length=10),
                   seq(11,20,length=10),
                   seq(21,30,length=10),
                   seq(31,40,length=10),
                   seq(41,50,length=10),
                   seq(51,60,length=10),
                   seq(61,70,length=10),
                   seq(71,80,length=10),
                   seq(81,90,length=10),
                   seq(91,100,length=10)),
          symm=F,
          symkey=F,
          symbreaks = T,
          cexCol = 1,
          col=colorRampPalette(c("white",
                                  "forestgreen"))(100))
  }
  dev.off()
}
