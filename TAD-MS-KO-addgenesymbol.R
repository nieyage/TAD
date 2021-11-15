rm(list = ls())
setwd("D:/project/TAD")

count1<-read.table("TAD-rawcounts.txt")
count2<-read.table("star_counts.txt")
count2<-count2[,c(1:3,7:9,13:18)]
####manager matrix#####
ENSEMBL <- gsub("\\_\\d*", "", count1$gene)
ENSEMBL <- gsub("\\.\\d*", "", ENSEMBL)
rownames(count1)<-ENSEMBL

library('biomaRt')
library("curl")
library(ggrepel)
library(dplyr)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))#####Ð¡Êó£ºmmusculus_gene_ensembl
my_ensembl_gene_id<-row.names(count1)
options(timeout = 4000000)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values= my_ensembl_gene_id, mart = mart)
head(mms_symbols)
head(count1)
rownames(mms_symbols)<-mms_symbols$ensembl_gene_id
mms_symbols<-mms_symbols[rownames(count1),]
head(mms_symbols)
count1$gene<-mms_symbols$external_gene_name
head(count1)
rownames(count1)<-count1$gene
###count2
ENSEMBL <- gsub("\\_\\d*", "", rownames(count2))
ENSEMBL <- gsub("\\.\\d*", "", ENSEMBL)
rownames(count2)<-ENSEMBL
head(count2)
my_ensembl_gene_id<-row.names(count2)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values= my_ensembl_gene_id, mart = mart)
head(mms_symbols)
head(count2)
rownames(mms_symbols)<-mms_symbols$ensembl_gene_id
mms_symbols<-mms_symbols[rownames(count2),]
head(mms_symbols)
count2$gene<-mms_symbols$external_gene_name
head(count2)
count1<-count1[!duplicated(count1$gene),]
count2<-count2[!duplicated(count2$gene),]
gene<-intersect(count1$gene,count2$gene)
gene<-na.omit(gene)
count1<-count1[which(count1$gene%in%gene),]
count2<-count2[which(count2$gene%in%gene),]

count<-merge(count1,count2,by="gene")
str(count)
count<-count[-1,]
tail(count)
rownames(count)<-count$gene
head(count)
count<-count[,c(14:16,2:4,17:19,8:10,20:22,11:13,23:25,5:7 )]
colnames(count)<-c("H1_Day0_rep1","H1_Day0_rep2","H1_Day0_rep3","MS_Day0_rep1","MS_Day0_rep2","MS_Day0_rep3",
                   "H1_Day2_rep1","H1_Day2_rep2","H1_Day2_rep3","MS_Day2_rep1","MS_Day2_rep2","MS_Day2_rep3",
                   "H1_Day5_rep1","H1_Day5_rep2","H1_Day5_rep3","MS_Day5_rep1","MS_Day5_rep2","MS_Day5_rep3",
                   "H1_Day10_rep1","H1_Day10_rep2","H1_Day10_rep3","MS_Day10_rep1","MS_Day10_rep2","MS_Day10_rep3")
write.table(count,"H1_MS_genesymbol_rawcount.txt")
