#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)
log <- file(args[5], open="wt")
sink(log, type="message")

id <- args[1]
annot1 <- args[2]
annot2 <- args[3]
output <- args[4]

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))


id<-fread(id)
a1<-fread(annot1)
a2<-fread(annot2)

id1<-merge(id,a1,by.x='V8',by.y='V1',all.x=T)
id2<-merge(id,a2,by.x='V8',by.y='V1',all.x=T)

id11 <- group_by(id1, V8) %>%
  summarise_all(funs(paste(unique(.), collapse = "|")))
id22 <- group_by(id2, V8) %>%
  summarise_all(funs(paste(unique(.), collapse = "|")))
  
res<-merge(id11,id22[,c(1,9,10)],by.x='V8',by.y='V8')
res$V8<-NULL
colnames(res) <- c('chr_bin1','start_bin1','end_bin1',
'chr_bin2','start_bin2','end_bin2','count',
'gene_bin1','exon_bin1','gene_bin2','exon_bin2')

fwrite(res,file=output, row.names=F, col.names=T,quote=F, sep='\t')
