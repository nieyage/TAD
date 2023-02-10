library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")

rm(list = ls())
setwd("D:/project/TAD/RNAseq/rawcounts")
H1<-read.table("H1_batch3.txt")
H1<-H1[,c(1:3,7:9,13:18)]
SA<-read.table("SA-KO-batch2.txt",header=T)
MS<-read.table("MS-KO-batch1.txt",header=T)

#modified rownames#
rownames(H1) <- gsub("\\.\\d*", "", rownames(H1))
MS<-MS[-1:-5,]
SA<-SA[-1:-5,]
MS$gene <- gsub("\\_\\d*", "", MS$gene)
MS$gene <- gsub("\\.\\d*", "", MS$gene)

SA$gene <- gsub("\\_\\d*", "", SA$gene)
SA$gene <- gsub("\\.\\d*", "", SA$gene)
head(SA)
# merge three matrix#
counts<-merge(SA,MS,by="gene")
H1$gene<-rownames(H1)
counts<-merge(counts,H1,by="gene")
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
rownames(counts)<-counts$gene
counts<-counts[mms_symbols$ensembl_gene_id,]
counts$genesymbol<-mms_symbols$external_gene_name
head(counts)
counts<-counts[!duplicated(counts$genesymbol),]
rownames(counts)<-counts$genesymbol
counts<-counts[,c(2:39)]
counts<-counts[,c(27:38,15:26,1:14)]
counts<-counts[,c(1:15,19:24,16:18,25:29,33:38,30:32)]
colnames(counts)<-c("H1_Day0_rep1_batch3","H1_Day0_rep2_batch3","H1_Day0_rep3_batch3",
                    "H1_Day2_rep1_batch3","H1_Day2_rep2_batch3","H1_Day2_rep3_batch3",
                    "H1_Day5_rep1_batch3","H1_Day5_rep2_batch3","H1_Day5_rep3_batch3",
                    "H1_Day10_rep1_batch3","H1_Day10_rep2_batch3","H1_Day10_rep3_batch3",
                    "MS_day0_rep1_batch1","MS_day0_rep2_batch1","MS_day0_rep3_batch1",
                    "MS_day2_rep1_batch1","MS_day2_rep2_batch1","MS_day2_rep3_batch1",
                    "MS_day5_rep1_batch1","MS_day5_rep2_batch1","MS_day5_rep3_batch1",
                    "MS_day10_rep1_batch1","MS_day10_rep2_batch1","MS_day10_rep3_batch1",
                    "H1_Day5_rep1_batch2","MS_Day5_rep1_batch2","SA_Day0_rep1_batch2",
                    "SA_Day0_rep2_batch2","SA_Day0_rep3_batch2","SA_Day2_rep1_batch2",
                    "SA_Day2_rep2_batch2","SA_Day2_rep3_batch2","SA_Day5_rep1_batch2",
                    "SA_Day5_rep2_batch2","SA_Day5_rep3_batch2","SA_Day10_rep1_batch2",
                    "SA_Day10_rep2_batch2","SA_Day10_rep3_batch2")
write.csv(counts,"TAD_all_rawcounts.csv")
write.table(counts,"TAD_all_rawcounts.txt")
