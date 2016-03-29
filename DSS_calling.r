# not quite working automatically, need to enter manually
# DSS DMR calling, make sure files are converted into right format using DSS_file_prep
options(echo=T)
args=commandArgs(trailingOnly=T)
print(args)

# install DSS
# source("http://bioconductor.org/biocLite.R")
# biocLite("DSS")

# use these arguments in order when running on server
context = args[1]
pvalue = args[2]
group1 = args[3]
group2 = args[4]

library(DSS)

#read in files, edited from DSS_file_prep.sh.
files <- dir(pattern = paste0(context,".output"))

# setup treatment lists
control <- list(files[5:8])
el1hr <- list(files[9:12])
el1hr_rec <- list(files[1:4])

# define group 1 and group 2 for comparison
group1 = control
group2 = el1hr

# read input files in DSS format (chr, pos, N, X)
dat1.1 <- read.delim(unlist(group1)[1])
dat1.2 <- read.delim(unlist(group1)[2])
dat1.3 <- read.delim(unlist(group1)[3])
dat1.4 <- read.delim(unlist(group1)[4])

dat2.1 <- read.delim(unlist(group2)[1])
dat2.2 <- read.delim(unlist(group2)[2])
dat2.3 <- read.delim(unlist(group2)[3])
dat2.4 <- read.delim(unlist(group2)[4])

# setup bsseq object
BSobj <- makeBSseqData(list(dat1.1,dat1.2,dat1.3,dat1.4,
			    dat2.1,dat2.2,dat2.3,dat2.4), 
	      sampleNames=c("C1","C2","C3","C4",
		 	    "N1","N2","N3","N4"))

BSobj <- makeBSseqData(list(dat1.1,dat1.2,dat1.3,
                            dat2.1,dat2.2,dat2.3),
              sampleNames=c("C1","C2","C3",
                            "N1","N2","N3"))

BSobj <- makeBSseqData(list(dat1.1,dat1.2,dat2.1,dat2.2),
              sampleNames=c("C1","C2","N1","N2"))

##################################
# DML testing without smoothing
dmlTest <- DMLtest(BSobj, group1=c("C1","C2","C3","C4"), group2=c("N1","N2","N3","N4"), smoothing=TRUE)
####################################

#############################################################
# perform dml testing with smoothing
dmlTest <- DMLtest(BSobj, group1=c("C1","C2","C3"), group2=c("N1","N2","N3"), smoothing=TRUE)
###############################################################
dmlTest <- DMLtest(BSobj, group1=c("C1","C2"), group2=c("N1","N2"), smoothing=TRUE)


# identify DMRs based on dmltesting
dmrs <- callDMR(dmlTest, delta=0.1, minlen=50, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=pvalue)
dmrs <- callDMR(dmlTest, delta=0.1, minlen=50, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=0.05)

# write output. as bedfile so can use bedtools for later analyses
file1=paste0(group1, "vs", group2, "_",context, "_", "p=", pvalue, ".bed")

write.table(dmrs,file=file1,quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
