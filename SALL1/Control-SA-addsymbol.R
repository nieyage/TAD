library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")

rm(list = ls())
setwd("D:/project/TAD/RNAseq/rawcounts")

SA<-read.table("SA-KO-batch2.txt",header=T)
control<-read.table("Control-batch4.txt",header=T)
counts<-merge(control,SA,by="gene")
counts$gene <- gsub("\\_\\d*", "", counts$gene)
counts$gene <- gsub("\\.\\d*", "", counts$gene)
counts<-counts[-1:-5,]
str(counts)
#trans ID to genesymbol#
library('biomaRt')
library("curl")
library(ggrepel)
library(dplyr)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))# mmusculus_gene_ensembl
my_ensembl_gene_id<-counts$gene
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values= my_ensembl_gene_id, mart = mart)
head(mms_symbols)
head(counts$gene)
colnames(counts)[1]<-"ensembl_gene_id"
counts<-merge(counts,mms_symbols,by="ensembl_gene_id")
head(counts)
counts<-counts[!duplicated(counts$external_gene_name),]
rownames(counts)<-counts$external_gene_name
counts<-counts[,c(2:28)]
counts<-counts[,-14:-15]
counts<-counts[,c(1:3,7:12,4:6,13:16,20:25,17:19)]
colnames(counts)<-c("Control_Day0_rep1_batch4", "Control_Day0_rep2_batch4", "Control_Day0_rep3_batch4",
                    "Control_Day2_rep1_batch4", "Control_Day2_rep2_batch4", "Control_Day2_rep3_batch4",
                    "Control_Day5_rep1_batch4", "Control_Day5_rep2_batch4", "Control_Day5_rep3_batch4",
                    "Control_Day10_rep1_batch4","Control_Day10_rep2_batch4","Control_Day10_rep3_batch4",
                    "SA_Day0_rep4_batch4",
                    "SA_Day0_rep1_batch2", "SA_Day0_rep2_batch2", "SA_Day0_rep3_batch2",
                    "SA_Day2_rep1_batch2", "SA_Day2_rep2_batch2", "SA_Day2_rep3_batch2",
                    "SA_Day5_rep1_batch2", "SA_Day5_rep2_batch2", "SA_Day5_rep3_batch2",
                    "SA_Day10_rep1_batch2","SA_Day10_rep2_batch2","SA_Day10_rep3_batch2")

write.csv(counts,  "TAD_SA_control_rawcounts.csv")
write.table(counts,"TAD_SA_control_rawcounts.txt")
