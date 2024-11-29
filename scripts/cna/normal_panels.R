#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)
log <- file(args[6], open="wt")
sink(log, type="message")

suppressPackageStartupMessages(library(ichorCNA))

createPanelOfNormals(filelist=args[1],
                     gcWig=args[2],
                     mapWig=args[3],
                     centromere=args[4],
                     outfile=args[5],
                     genomeBuild = "hg38")