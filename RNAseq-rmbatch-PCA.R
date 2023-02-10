library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("D:/project/TAD/RNAseq/")
counts<-read.table("./rawcounts/TAD_all_rawcounts.txt")
head(counts)
# create pesudo replicates for batch2#
counts$H1_Day5_rep2_batch2<-counts$H1_Day5_rep1_batch2
counts$MS_Day5_rep2_batch2<-counts$MS_Day5_rep1_batch2
str(counts)
counts<-counts[,c(1:24,25,39,26,40,27:38)]
#remove batch effect#
condition <- factor(c(rep("H1",12),rep("MS",12),rep("H1",2),rep("MS",2),rep("SA",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day5",4),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))
batch <- factor(c(rep("batch3",12),rep("batch1",12),rep("batch2",16)))

colData <- data.frame(row.names=colnames(counts),
                      condition=condition,timepoint=timepoint,
                      batch=batch)
colData
colData$condition <- relevel(colData$condition, ref="H1")

countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData<-na.omit(countData)
dds<-DESeqDataSetFromMatrix(countData,colData, 
                            formula(~timepoint+condition+batch)) 
dds <- estimateSizeFactors(dds)
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
#rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
#pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
#hc <- hclust(t(rlogMat) )

#PCA not remove batch
library(ggplot2)
vsd <- vst(dds)
head(vsd)
pca_data <- plotPCA(vsd, intgroup=c("batch","timepoint","condition"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p<-ggplot(pca_data, aes(PC1, PC2, shape =condition,color=timepoint)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
geom_text(aes(label=name),vjust=2)
p+theme_bw()
## PCA remove batch 
library(limma)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), c(colData$batch))
pca <- plotPCA(vsd, intgroup=c("batch","timepoint","condition"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca, "percentVar"))
p<-ggplot(pca[1:24,], aes(PC1, PC2, shape =condition,color=timepoint)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
geom_text(aes(label=name),vjust=2)
p
p+theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

 head(assay(vsd))
pcacount<-assay(vsd)
head(pcacount)
write.table(pcacount,"RNAseq-PCAcount.csv")
