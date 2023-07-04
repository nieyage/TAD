library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("D:/project/TAD/RNAseq/")
counts<-read.table("./rawcounts/TAD_SA_control_rawcounts.txt")
head(counts)
# create pesudo replicates for batch2#
counts$SA_Day0_rep5_batch4<-counts$SA_Day0_rep4_batch4
str(counts)
counts<-counts[,c(1:13,26,14:25)]
#remove batch effect#
condition <- factor(c(rep("Control",12),rep("SA",14)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",5),rep("Day2",3),rep("Day5",3),rep("Day10",3)))
batch <- factor(c(rep("batch4",14),rep("batch2",12)))
colData <- data.frame(row.names=colnames(counts),
                      condition=condition,timepoint=timepoint,
                      batch=batch)
colData
colData$condition <- relevel(colData$condition, ref="Control")

countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData<-na.omit(countData)
dds<-DESeqDataSetFromMatrix(countData,colData, 
                            formula(~timepoint+condition+batch+condition:timepoint)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dds <- DESeq(dds)
res <- results(dds)
#condition DEG#
SA_Control_Day2<-results(dds, contrast=list( c("condition_SA_vs_Control",
                                               "timepointDay2.conditionSA") ))
SA_Control_Day2_diff_gene <-rownames(subset(SA_Control_Day2, 
                                            padj < 0.05 & abs(log2FoldChange)> 1))
#Day2 vs Day0
Day2_0 <-rownames(subset(results(dds,contrast = c("timepoint","Day2","Day0")),
                         padj < 0.05 & abs(log2FoldChange)> 1)) 

#plot all Day2_Day0 DEG:

targetcount<-normalized_counts[Day2_0,c(1:6,18:20)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
condition <- factor(c(rep("Control",6),rep("SALL1 KD",3)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",6)))

annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)
rownames(annotation_col) = factor(colnames(count))

ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3"),
  Group = c("Control"="#F39B7FB2","SALL1 KD" ="#8491B4B2"))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
#All Day0_2 gene#
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         cutree_rows =5,
         clustering_method ='ward.D2',
         cutree_cols = 2,
         #annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=F,show_colnames=F)
###change signif gene heatmap 
SA_gene <-rownames(subset(SA_Control_Day2,padj < 0.05 & log2FoldChange > 1))
Control_gene <-rownames(subset(SA_Control_Day2,padj < 0.05 & log2FoldChange < -1))

#Day2 vs Day0
Day2 <-rownames(subset(results(dds,contrast = c("timepoint","Day2","Day0")),padj < 0.05 & log2FoldChange> 1)) 
Day0 <-rownames(subset(results(dds,contrast = c("timepoint","Day2","Day0")),padj < 0.05 & log2FoldChange< -1)) 

cluster1<-intersect(Day2,Control_gene)
cluster2<-intersect(Day0,SA_gene)
write.csv(cluster1,"SALL1-Day2vsDay0-down-up-down.csv")
write.csv(cluster2,"SALL1-Day2vsDay0-up-down-up.csv")

gene<-c(cluster1,cluster2)
targetcount<-normalized_counts[gene,c(1:6,18:20)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
#All Day0_2 gene#
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         cutree_rows =2,
         clustering_method ='ward.D2',
         cutree_cols = 2,
         #annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=F,show_colnames=F)

####GO and KEGG for down-up-down gene:
library(clusterProfiler)
library(org.Hs.eg.db)
down-up-down gene:
  pdf("down-up-down-gene-GO.pdf")
gene.df <- bitr(cluster1, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"down-up-down-gene-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"down-up-down-gene-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"down-up-down-gene-CC.csv")
dev.off()
######KEGG##########
pdf("down-up-down-gene_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"down-up-down-gene-KEGG.csv")

dev.off()

####GO and KEGG for up-down-up gene:

pdf("up-down-up-gene-GO.pdf")
gene.df <- bitr(cluster2, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"up-down-up-gene-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"up-down-up-gene-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"up-down-up-gene-CC.csv")
dev.off()
######KEGG##########
pdf("up-down-up-gene_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"up-down-up-gene-KEGG.csv")

dev.off()


