#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)


data_list <- args[1]
data_merged <- args[2]

data_list<-as.vector(read.table(data_list,header=F))$V1

data1<-data.frame()
for (i in data_list) {
tmp<-read.delim(i)
tmp$Sample<-gsub('.*/','',gsub('-plasma-tumor-shared.cna.txt','',i))
data1<-rbind(data1,tmp)
}


res <- data1[with(data1, order(data1$Sample, data1$Tumor.LogR,decreasing=T)), ]
res<-res[,c(ncol(res),1:(ncol(res)-1))]
write.table(res, file = data_merged, col.names=T, row.names=FALSE, sep="\t",quote=F)