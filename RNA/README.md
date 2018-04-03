# RNA-seq scripts repository

#### RNAseq_bam_to_100bpwigs.sh
Produce 100bp windows (100bp.bed) of RNAseq coverage from BAMs across annotations of interest.

#### RNAseq_bam_to_bedgraph.sh	
Produce coverage data from BAM files for RNAseq or ChIP data in bedgraph format, subsequently converting to bigWig files (IGV browsing). 

#### RNAseq_featureCounts.sh
Perform feature counts on aligned reads on annotation file of interest for differential expression analysis.

#### RNAseq_v0.1.sh
Master script for performing quality trimming and subread alignment of RNA-seq reads. Sorted BAM as output for down-stream analysis.

#### edgeR.r
Template script to run edgeR for differential gene expression calling.

#### rel_expression_plots.r
Get raw read coverage across features of interest from RNA WIG files (see RNAseq_bam_to_100bpwigs).

