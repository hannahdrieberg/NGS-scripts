#!/bin/bash
set -eu

# convert bismark cov files to IGV format, then produce TDF files of that data.
# Use in directory with all cov files of interest.
# Can produce TDF for depth or prop. methylation

for FILE in *bismark.cov
do
cat -n $FILE | awk -v OFS="\t" '{print $2, $3-1, $4, $1, $5}' > ${FILE%%.bismark.cov}.igv
java -Xmx26g -Djava.awt.headless=true -jar /home/diep/bin/IGVTools/igvtools.jar toTDF ${FILE%%.bismark.cov}.igv ${FILE%%.bismark.cov}.tdf /home/diep/Araport11/Araport11.genome
done

# Depth of coverage
#cov=read.delim(flist[i],head=F)
#cov$V2=cov$V2-1
#cov$id=seq(1:nrow(cov))
#cov$V7=cov$V5+cov$V6
#cov=cov[,c(1,2,3,7,8)]
