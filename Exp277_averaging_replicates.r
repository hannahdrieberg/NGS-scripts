# R script for Exp277 Control vs Drought Samples
# Mapped DMR position to genomic features = loci that are near/contain DMRs (any context)
# Extracted methylation information within these loci, and upstream, and downstream; using wig files.
# Read this data in and average replicates together to create a dataframe to start visualising
# Extracting sample information based on file naming - THIS IS SUBJECT TO CHANGE EVERY EXPERIMENT
# Set working directory to where files are 
# make sure to have files of a particular feature there or make sure naming is appropriate
a <- dir(pattern = "_gene")
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
# substr replicate info BASED ON SAMPLE NAME
x <- as.numeric(regexec(pattern = "7-", text = a[i])) 
y <- as.numeric(regexec(pattern = "_C", text = a[i]))
z <- as.numeric(regexec(pattern = "_D", text = a[i]))
input$sample <- substr(x = as.character(a[i]), start = x+2, stop = y-1) 
input$context <- substr(x = as.character(a[i]), start = y+1, stop = z-1)
out <- rbind(input, out)
}
out$Treatment <- factor(ifelse(test = out$sample == "D7" | 
                          out$sample == "D8" |
                          out$sample == "D9", "Drought", "Control"))


# average reps
out.reps <- ddply(out, .(V1,V2,V3,V9,V10,V11,V12,V13,V14,V15,real.dist,Treatment,context), summarise,
             Prop_met = mean(V4),
             Coverage = mean(V7),
             C_count = mean(V8),
             n = sum(!is.na(V4)) # a few funny things happening
             )
# Gene Body
out.body <- subset(out.reps,out.reps$real.dist == 0)
rel.dist=matrix(ifelse(out.body[,9]=="-",((out.body[,6]-(out.body[,2]))/(out.body[,6]-out.body[,5]))*1000,(((out.body[,2])-out.body[,5])/(out.body[,6]-out.body[,5]))*1000),ncol=1)
out.body=cbind(out.body,rel.dist)
fixy=ifelse(out.body$rel.dist < 0,0,ifelse(out.body$rel.dist>1000,1000,out.body$rel.dist))
fixy <- round(fixy, digits = -2)
out.body$rel.dist=fixy
out.v1<-out.body
out.v2=dcast(out.body, V12 ~ factor(rel.dist)+Treatment+context, mean, value.var = "Prop_met")