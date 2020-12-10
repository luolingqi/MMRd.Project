#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
   stop("Exactly 1 argument must be supplied: project", call.=FALSE)
   }

library("SigProfilerPlottingR")
project = args[1]

plot_path <- paste0("/data/",project,"/output/plot/")
if (! file.exists(plot_path)){
   dir.create(plot_path)
   }

out_dir <- paste0("/data/",project,"/output")

DBS_type=c("1248","186","2976","78")
for (i in DBS_type) {
    plotDBS(paste0(out_dir,"/DBS/",project,".DBS",i,".all"), file.path(out_dir,"plot/"), project, i, percentage=TRUE)
    }

SBS_type=c("1536","384","6144","24","6","96")
for (i in SBS_type) {
    plotSBS(paste0(out_dir,"/SBS/",project,".SBS",i,".all"), file.path(out_dir,"plot/"), project, i, percentage=TRUE)
    }

ID_type=c("28","415","83","8628","96")
for (i in ID_type) {
    plotID(paste0(out_dir,"/ID/",project,".ID",i,".all"), file.path(out_dir,"plot/"), project, i, percentage=TRUE)
    }