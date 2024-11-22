#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)

suppressPackageStartupMessages(library(ichorCNA))

subclonal= args[6]

if (args[6]=='plasma') {
estimateSc=F
} else {
estimateSc=T
}

run_ichorCNA(id=args[1],
             tumor_wig=args[2],
             ploidy=args[3],
             normal=args[4],
             maxCN=as.numeric(args[5]),
             estimateScPrevalence=estimateSc,
             scStates=args[7],
             normal_panel=args[8],
             outDir=args[9],
             cores=args[10],
             gcWig=args[11],
             mapWig=args[12],
             centromere = args[13],
             minMapScore = args[14],
             chrs="c(1:22)",
             chrTrain= "c(1:22)",
             fracReadsInChrYForMale = 0.002,
             genomeBuild = "hg38")