#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)

bed <- args[1]
tmp_dir <- args[2]
cores <- as.numeric(args[3])

setwd(tmp_dir)
suppressPackageStartupMessages(rgumbo::splitBed(bed,chunks=cores,writeBed = T,verbose = F))