#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rgumbo))
suppressPackageStartupMessages(library(tidyverse))

fb_table <- args[1]
vs_table <- args[2]

#LoFreq
lf_tumor <- args[3]
lf_plasma <- args[4]
#Freebayes
fb_tumor <- args[5]
fb_plasma <- args[6]
#Varscan2
vs_tumor <- args[7]
vs_plasma <- args[8]
#Outputs
call<-args[9]
sum_stats<-args[10]

fb<-fread(fb_table)
vs<-fread(vs_table)

if (nrow(vs)>0 & nrow(fb)>0) {

vs$id<-paste(vs$Chrom,vs$Pos,sep=':')
fb$id<-paste(fb$Chrom,fb$Pos,sep=':')
m1<-merge(vs,fb,by.x='id',by.y='id',all.x=T,all.y=T,sort=F)

colnames(m1)[c(5:7,11:13)]<-c('Varscan2_control','Varscan2_plasma','Varscan2_tumor','Freebayes_control','Freebayes_plasma','Freebayes_tumor')
m1$Calling_score<-ifelse(apply(m1[,c(4:6,10:12)],1,function(x) sum(is.na(x)))==0,2,1)
m2<- m1 %>% mutate(Chr = coalesce(Chrom.x ,Chrom.y)) %>%
mutate(Pos = coalesce(Pos.x ,Pos.y))  %>%
mutate(Control = coalesce(Varscan2_control ,Freebayes_control))  %>%
mutate(Plasma = coalesce(Varscan2_plasma ,Freebayes_plasma))  %>%
mutate(Tumor = coalesce(Varscan2_tumor ,Freebayes_tumor))  %>%
select(id,Chr, Pos, Control, Plasma, Tumor, Calling_score)
m3<-m2[order(m2$Calling_score,decreasing=T),]

vcf=lf_tumor
vcf<- fread(vcf,skip=grep('#CHROM',readLines(vcf))-1)
colnames(vcf)[1]<-'Chr'
vcf$id<-paste(vcf$Chr,vcf$POS,sep=':')
m4<-merge(m3,vcf[,c('REF','ALT','id')],by.='id',by.y='id',all.x=T,all.y=T,sort=F)
m4$lf_tumor<-ifelse(!is.na(m4$REF),paste(m4$REF,m4$ALT,sep='/'),NA)
m4$REF<-NULL;m4$ALT<-NULL

vcf=lf_plasma
vcf<- fread(vcf,skip=grep('#CHROM',readLines(vcf))-1)
colnames(vcf)[1]<-'Chr'
vcf$id<-paste(vcf$Chr,vcf$POS,sep=':')
m5<-merge(m4,vcf[,c('REF','ALT','id')],by.='id',by.y='id',all.x=T,all.y=T,sort=F)
m5$lf_plasma<-ifelse(!is.na(m5$REF),paste(m5$REF,m5$ALT,sep='/'),NA)
m5$REF<-NULL;m5$ALT<-NULL

m5$Calling_score<-ifelse(apply(m5[,c(8:9)],1,function(x) sum(is.na(x)))<2,m5$Calling_score+1,m5$Calling_score)
m5$Calling_score<-ifelse(is.na(m5$Calling_score),1,m5$Calling_score)
m5$Control<-ifelse(is.na(m5$Pos),gsub('/.*','',m5$lf_tumor),m5$Control)
m5$Control<-ifelse(is.na(m5$Pos) & is.na(m5$Control),gsub('/.*','',m5$lf_plasma),m5$Control)
m5$Plasma <- ifelse(is.na(m5$Plasma),m5$lf_plasma,m5$Plasma)
m5$Plasma <- ifelse(is.na(m5$Plasma),0,m5$Plasma)
m5$Tumor <- ifelse(is.na(m5$Tumor),m5$lf_tumor,m5$Tumor)
m5$Tumor <- ifelse(is.na(m5$Tumor),0,m5$Tumor)

m6<-m5[order(m5$Calling_score,decreasing=T),c(1,4:7)]
m7<-m6[!grepl('chrUn|random',m6$id),]
}

fb_tumor<-SnpEff(fb_tumor,modifier=T)
fb_tumor <- group_by(fb_tumor, Genomic_coordinates) %>%
  summarise_all(funs(paste(unique(.), collapse = "|")))

fb_plasma<-SnpEff(fb_plasma,modifier=T)
fb_plasma <- group_by(fb_plasma, Genomic_coordinates) %>%
  summarise_all(funs(paste(unique(.), collapse = "|")))

vs_tumor<-SnpEff(vs_tumor,modifier=T)
vs_tumor <- group_by(vs_tumor, Genomic_coordinates) %>%
  summarise_all(funs(paste(unique(.), collapse = "|")))

vs_plasma<-SnpEff(vs_plasma,modifier=T)
vs_plasma <- group_by(vs_plasma, Genomic_coordinates) %>%
  summarise_all(funs(paste(unique(.), collapse = "|")))

annot<-rbind(fb_tumor[,c(1:8,10:11)],fb_plasma[,c(1:8,10:11)],vs_tumor[,c(1:8,10:11)],vs_plasma[,c(1:8,10:11)])
annot1<-annot[!duplicated(annot$Genomic_coordinates),]

res<-merge(m7,annot1,by.x='id',by.y='Genomic_coordinates',all.x=T)
res <- res[with(res, order(res$impact,res$id)), ]
res$chr<-gsub(':.*','',res$id)
res$pos<-gsub('.*:','',res$id)
res$id<-NULL
res<-res[,c(14,15,12:13,5:11,1:4)]
sum1<-data.frame(N.somatic.SNVs=c(sum(res$Plasma!=0),sum(res$Tumor!=0)),
                 N.shared.somatic.SNVs=c(sum(res$Plasma!=0 & res$Tumor!=0),''),
                 Perc.shared.SNVs=c(round(sum(res$Plasma!=0 & res$Tumor!=0)/sum(res$Plasma!=0),digits=2),''))
fwrite(res,file=call, row.names=F, col.names=T,quote=F, sep='\t')
fwrite(sum1,file=sum_stats, row.names=F, col.names=T,quote=F, sep='\t')