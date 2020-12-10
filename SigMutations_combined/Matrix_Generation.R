#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
   stop("Exactly 3 arguments must be supplied: 1) Project_name; 2) reference_name; 3)path_to_input_vcf.n", call.=FALSE)
   }

library("SigProfilerMatrixGeneratorR")
project = args[1]
reference = args[2]
path_vcf = args[3]

assign(paste0("matrices_",project),SigProfilerMatrixGeneratorR(project, reference, path_vcf, plot=F, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100))


#matrices_DTH_gain_exome <- SigProfilerMatrixGeneratorR("DTH_gain_exome", "mm10", "/data/P_vs_DTH_gain_exome/", plot=F, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
