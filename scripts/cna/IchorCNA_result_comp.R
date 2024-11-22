#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)
log <- file(args[5], open="wt")
sink(log, type="message")

plasma <- args[1]
tumor <- args[2]
annot <- args[3]
output <- args[4]

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

ann<-fread(annot)
tmp <- GRanges(ann$V1, ranges =IRanges(ann$V2,ann$V3))

#Plasma seg
cna_p<-as.data.frame(fread(plasma))
cna_p <- cna_p[,c(1:3,6,8:10)]
cna_p$chr<-paste0('chr',cna_p$chr)
g1<-GRanges(cna_p$chr, ranges =IRanges(cna_p$start,cna_p$end))
type3 = findOverlaps(query = tmp, subject = g1)
cap = data.frame(cna_p[subjectHits(type3),], ann[queryHits(type3),])
cap<-cap[,c(1:7,11)]
colnames(cap)<-c('Chr','Start','End','Plasma.LogR','Plasma.Copy_Number','Plasma.Call','Plasma.LogR_Copy_Number','Gene_name')
cap1<-cap[!is.na(cap$Plasma.LogR_Copy_Number),]
cap1$pos<- paste(cap1$Chr,cap1$Start,cap1$End,sep='-')
cap2<- as.data.frame(group_by(cap1, pos) %>%
  summarise_all(funs(paste(unique(.), collapse = "|"))))
cap2$Plasma.Call <- ifelse(cap2$Plasma.Call=='NEUT','Copy Neutral', ifelse(cap2$Plasma.Call=='HETD','Hemizygous Deletions',
ifelse(cap2$Plasma.Call=='GAIN','Copy Gain',ifelse(cap2$Plasma.Call=='HLAMP','High-level Amplification',
ifelse(cap2$Plasma.Call=='AMP','Amplification',ifelse(cap2$Plasma.Call=='HOMD','Homozygous Deletions State',cap2$Plasma.Call))))))

#Tumor seg
cat<-as.data.frame(fread(tumor))
cat <- cat[,c(1:3,6,8:10)]
cat$chr<-paste0('chr',cat$chr)
colnames(cat)<-c('Chr','Start','End','Tumor.LogR','Tumor.Copy_Number','Tumor.Call','Tumor.LogR_Copy_Number')
cat1<-cat[!is.na(cat$Tumor.LogR_Copy_Number),]
cat1$pos<- paste(cat1$Chr,cat1$Start,cat1$End,sep='-')
cat2<- as.data.frame(group_by(cat1, pos) %>%
  summarise_all(funs(paste(unique(.), collapse = "|"))))
cat2$Tumor.Call <- ifelse(cat2$Tumor.Call=='NEUT','Copy Neutral', ifelse(cat2$Tumor.Call=='HETD','Hemizygous Deletions',
ifelse(cat2$Tumor.Call=='GAIN','Copy Gain',ifelse(cat2$Tumor.Call=='HLAMP','High-level Amplification',
ifelse(cat2$Tumor.Call=='AMP','Amplification',ifelse(cat2$Tumor.Call=='HOMD','Homozygous Deletions State',cat2$Tumor.Call))))))

#Comparison
ca<-merge(cat2,cap2,by.x='pos',by.y='pos')
ca$pos<-NULL
ca<-ca[,c(1:3,15,4:7,11:14)]
colnames(ca)[1:3]<-c('Chr','Start','End')
fwrite(ca,file=output, row.names=F, col.names=T,quote=F, sep='\t')
