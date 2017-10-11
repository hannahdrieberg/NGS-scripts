# NGS-scripts
Repository for scripts used for analysing next generation sequencing data (described in order below). 

### 100bp_dmr_merge.r
Supplemental script for *100bp_dmrs.v0.1.sh* for collapsing/summarising DMRs

### 100bp_dmrs.v0.1.sh
This script uses 100bp windowed methylation data to call differences across the genome in a defined context with required coverage and a defined difference. These are performed in a pairwise manner between all samples of interest.

### 100bp_heatmap.sh
Uses the output from *100bp_dmrs.v0.1.sh* to get mC from bed files (individual mC resolution).

### 100bp_wig_to_dmrs.r
First step of *100bp_dmrs.v0.1.sh* that makes pairwise comparisons of windows showing a defined difference in a given context.

### BS-SNPer.sh		
Script to perform SNP calling from aligned bisulfite converted reads. Use on sorted BAM file.

### C_context_window_SREedits.pl	
SRE perl script to create 100bp windows from bismark methylation extractor report files. 

### DSS_calling.r
Script for performing DSS DMR calling on re-formatted (using DSS_file_prep.r) bed.cov file output.

### DSS_file_prep.r
Script for re-formatting bed.cov file methylation output for input to DSS.

### RNAseq_bam_to_100bpwigs.sh
Produce 100bp windows (100bp.bed) of RNAseq coverage from BAMs across annotations of interest.

### RNAseq_bam_to_bedgraph.sh	
Produce coverage data from BAM files for RNAseq or ChIP data in bedgraph format, subsequently converting to bigWig files (IGV browsing). 

### RNAseq_featureCounts.sh
Perform feature counts on aligned reads on annotation file of interest for differential expression analysis.

### RNAseq_v0.1.sh
Master script for performing quality trimming and subread alignment of RNA-seq reads. Sorted BAM as output for down-stream analysis.

### TruSeq-adapters.fa
FASTA file containing illumina adapter sequences for scythe step in RNA-seq alignment.

### araport11_assemble.sh

### average_cov.sh
Calculate average depth using samtools depth on sorted BAM file.

### bed_to_rel_dist.sh

### chip-seq_v0.1.sh

### conversion_rate_check.sh

### dmr_merge.r

### gene_to_gene_anno.sh

### genome_wigs_heatmap.r

### merge_covs.r

### merge_wigs.r

### met-sign.sh

### methylation_tiling.sh

### rel_expression_plots.r

### rel_methylation_plots-v2.r

### rel_methylation_plots.r

### repeat_methylation_plots.r

### scatman_smooth.sh

### smooth_scat.r

### wgbs_cov_to_TDF.txt

### wgbs_pipelinev0.4.sh

### wgbs_wig_PCA.r

