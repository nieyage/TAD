library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("D:/project/TAD/RNAseq/MS-H1")
counts<-read.table("../rawcounts/TAD_all_rawcounts.txt")
head(counts)
# create pesudo replicates for batch2#
counts$H1_Day5_rep2_batch2<-counts$H1_Day5_rep1_batch2
counts$MS_Day5_rep2_batch2<-counts$MS_Day5_rep1_batch2
str(counts)
counts<-counts[,c(1:24,25,39,26,40,27:38)]
counts<-counts[,1:28]
#remove batch effect#
condition <- factor(c(rep("H1",12),rep("MS",12),rep("H1",2),rep("MS",2)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day5",4)))
batch <- factor(c(rep("batch3",12),rep("batch1",12),rep("batch2",4)))
group <- factor(paste0(condition, timepoint))

colData <- data.frame(row.names=colnames(counts),
                      condition=condition,timepoint=timepoint,
                      batch=batch,group=group)
colData
colData$condition <- relevel(colData$condition, ref="H1")

countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData<-na.omit(countData)
dds<-DESeqDataSetFromMatrix(countData,colData, 
                            formula(~timepoint+condition+batch+condition:timepoint)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
dds <- DESeq(dds)
res <- results(dds)
#condition DEG#
MS_H1_Day2<-results(dds, contrast=list( c("condition_MS_vs_H1",
                                          "timepointDay2.conditionMS") ))
MS_H1_Day2_diff_gene <-rownames(subset(MS_H1_Day2, 
                                       padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(MS_H1_Day2,"MS_H1_Day2.csv")
#volcano plot#
pdf("MS-vs-H1-Day2-volcano.pdf")
data<-as.data.frame(MS_H1_Day2)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'MS-KO','Control'),
                     'Stable')
table(data$change)
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-12,12)) +
  ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)",
       title=paste("Day2:MS-KO vs Control DEG","\n","Control:553 genes; MS-KO:513 genes;",sep="")
       )+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%MS_H1_Day2_diff_gene),]
targetcount<-targetcount[,c(4:6,16:18)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("Day2-MS-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()

###Day0 
MS_H1_Day0<-results(dds, 
                    contrast=list( c("condition_MS_vs_H1") ))
MS_H1_Day0_diff_gene <-rownames(subset(MS_H1_Day0, 
                                       padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(MS_H1_Day0,"MS_H1_Day0.csv")
#volcano plot#
pdf("MS-vs-H1-Day0-volcano.pdf")
data<-as.data.frame(MS_H1_Day0)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'MS-KO','Control'),
                     'Stable')
table(data$change)
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-12,12)) +
  ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)",
       title=paste("Day0:MS-KO vs Control DEG","\n","Control:342 genes; MS-KO:266 genes;",sep=""))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

targetcount<-normalized_counts[which(rownames(normalized_counts)%in%MS_H1_Day0_diff_gene),]
targetcount<-targetcount[,c(1:3,13:15)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("Day0-MS-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()

#Day5
MS_H1_Day5<-results(dds, contrast=list( c("condition_MS_vs_H1",
                                          "timepointDay5.conditionMS") ))
MS_H1_Day5_diff_gene <-rownames(subset(MS_H1_Day5, 
                                       padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(MS_H1_Day5,"MS_H1_Day5.csv")
#volcano plot#
pdf("MS-vs-H1-Day5-volcano.pdf")
data<-as.data.frame(MS_H1_Day5)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'MS-KO','Control'),
                     'Stable')
table(data$change)
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-12,12)) +
  ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)",
       title=paste("Day5:MS-KO vs Control DEG","\n","Control:1656 genes; MS-KO:1788 genes;",sep=""))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%MS_H1_Day5_diff_gene),]
targetcount<-targetcount[,c(7:9,16:18)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("Day5-MS-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()

#Day10:
MS_H1_Day10<-results(dds, contrast=list( c("condition_MS_vs_H1",
                                           "timepointDay10.conditionMS") ))
MS_H1_Day10_diff_gene <-rownames(subset(MS_H1_Day10, 
                                        padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(MS_H1_Day10,"MS_H1_Day10.csv")
#volcano plot#
pdf("MS-vs-H1-Day10-volcano.pdf")
data<-as.data.frame(MS_H1_Day10)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'MS-KO','Control'),
                     'Stable')
table(data$change)
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$padj), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-12,12)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (padj)",
       title=paste("Day2:MS-KO vs Control DEG","\n","Control:1852 genes; MS-KO:1840 genes;",sep=""))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%MS_H1_Day10_diff_gene),]
targetcount<-targetcount[,c(7:9,19:21)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("Day10-MS-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()


library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
####GO and KEGG for MS-KO-DEG-Day2:
MS_H1_Day5_up <-subset(MS_H1_Day5, padj < 0.05 & log2FoldChange >1) 
MS_H1_Day5_up<-rownames(MS_H1_Day5_up)
MS_H1_Day5_up
pdf("MS_H1_Day5_up-GO.pdf")

gene.df <- bitr(MS_H1_Day5_up, fromType = "SYMBOL",
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
write.csv(ego,"MS_H1_Day5_upgene-BP.csv")


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
write.csv(ego,"MS_H1_Day5_upgene-MF.csv")

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
write.csv(ego,"MS_H1_Day5_upgene-CC.csv")
dev.off()
######KEGG##########
pdf("MS_H1_Day5_upgene-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hMS',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MS_H1_Day5_upgene-KEGG.csv")

dev.off()
########downgene############
MS_H1_Day5_down <-subset(MS_H1_Day5, padj < 0.05 & log2FoldChange < -1) 
str(MS_H1_Day5_down)
MS_H1_Day5_down<-rownames(MS_H1_Day5_down)
MS_H1_Day5_down:174 genes
pdf("MS_H1_Day5_downgene-GO.pdf")
gene.df <- bitr(MS_H1_Day5_down, fromType = "SYMBOL",
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
write.csv(ego,"MS_H1_Day5_downgene-GO-BP.csv")


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
write.csv(ego,"MS_H1_Day5-downgene-GO-MF.csv")

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
write.csv(ego,"MS_H1_Day5-downgene-GO-CC.csv")
dev.off()
######KEGG##########
pdf("MS_H1_Day5-downgene-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hMS',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MS_H1_Day5-downgene-KEGG.csv")
dev.off()


Day0_2 <-rownames(subset(results(dds,contrast = c("timepoint","Day0","Day2")),padj < 0.05 & log2FoldChange> 1)) 
Day0_5 <-rownames(subset(results(dds,contrast = c("timepoint","Day0","Day5")),padj < 0.05 & log2FoldChange> 1)) 
Day0_10 <-rownames(subset(results(dds,contrast = c("timepoint","Day0","Day10")),padj < 0.05 & log2FoldChange> 1)) 
Day0<-intersect(Day0_2,Day0_5)
Day0<-intersect(Day0,Day0_10)

Day2_0 <-rownames(subset(results(dds,contrast = c("timepoint","Day2","Day0")),padj < 0.05 & log2FoldChange> 1)) 
Day2_5 <-rownames(subset(results(dds,contrast = c("timepoint","Day2","Day5")),padj < 0.05 & log2FoldChange> 1)) 
Day2_10 <-rownames(subset(results(dds,contrast = c("timepoint","Day2","Day10")),padj < 0.05 & log2FoldChange> 1)) 
Day2<-intersect(Day2_0,Day2_5)
Day2<-intersect(Day2,Day2_10)
Day5_0 <-rownames(subset(results(dds,contrast = c("timepoint","Day5","Day0")),padj < 0.05 & log2FoldChange> 1)) 
Day5_2 <-rownames(subset(results(dds,contrast = c("timepoint","Day5","Day2")),padj < 0.05 & log2FoldChange> 1)) 
Day5_10 <-rownames(subset(results(dds,contrast = c("timepoint","Day5","Day10")),padj < 0.05 & log2FoldChange> 1)) 
Day5<-intersect(Day5_0,Day5_2)
Day5<-intersect(Day5,Day5_10)

Day10_0 <-rownames(subset(results(dds,contrast = c("timepoint","Day10","Day0")),padj < 0.05 & log2FoldChange> 1)) 
Day10_2 <-rownames(subset(results(dds,contrast = c("timepoint","Day10","Day2")),padj < 0.05 & log2FoldChange> 1)) 
Day10_5 <-rownames(subset(results(dds,contrast = c("timepoint","Day10","Day5")),padj < 0.05 & log2FoldChange> 1)) 
Day10<-intersect(Day10_0,Day10_2)
Day10<-intersect(Day10,Day10_5)

time_DEG<-c(Day0,Day5,Day10,Day2)
time_DEG<-time_DEG[!duplicated(time_DEG)]


targetcount<-normalized_counts[which(rownames(normalized_counts)%in%time_DEG),]
targetcount<-targetcount[,c(1:12,13:24)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
count<-na.omit(count)
p<-pheatmap(count,cluster_cols = F,cluster_rows = T,
            cutree_rows =6, 
            color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
            #cellwidth = 10, cellheight = 10,
            show_rownames=F,show_colnames=T)
row_cluster <- cutree(p$tree_row, k=6)
head(row_cluster)
Day0<-names(row_cluster[row_cluster==2])
Day2<-names(row_cluster[row_cluster==4])
Day5<-names(row_cluster[row_cluster==3])
Day10<-names(row_cluster[row_cluster==1])


count=t(scale(t(targetcount[c(Day0,Day2,Day5,Day10),]),scale = T,center = T))

pheatmap(count,cluster_cols = F,cluster_rows = F,
         cutree_rows =1, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)

setwd("D:/project/TAD/RNAseq/MS-H1/time point")
####GO and KEGG for MS-KO-Day0:
pdf("Day0-GO.pdf")
gene.df <- bitr(Day0, fromType = "SYMBOL",
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
write.csv(ego,"Day0-BP.csv")


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
write.csv(ego,"Day0-MF.csv")

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
write.csv(ego,"Day0-CC.csv")
dev.off()
######KEGG##########
pdf("Day0_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hMS',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"Day0-KEGG.csv")

dev.off()

write.csv(normalized_counts,"MS-H1-normalized_counts.csv")

#cell Net analysis
data<-read.csv("../CellNetAnalysis/classification_results_2022-03-22.csv",row.names = 1)
str(data)
data<-data[,1:24]
condition <- factor(c(rep("WT",12),rep("MEIS1 KO",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))

annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)
rownames(annotation_col) = factor(colnames(data))
pdf("MS-H1-CellNet.pdf",width = 10,height = 6)
pheatmap(data[c(1:12,14:15),],cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("black", "#48AF3A", "#E7E639"))(100),
         cellwidth = 15, cellheight = 10,
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         show_rownames=T,show_colnames=T)
dev.off()


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
        legend.title = element_blank())#+
  #geom_text(aes(label=name),vjust=2)
p
pdf("MS-H1-PCA.pdf",width = 8,height = 5)
p+theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()
head(assay(vsd))
pcacount<-assay(vsd)
head(pcacount)
write.table(pcacount,"MS-H1-RNAseq-PCAcount.csv")

#Stage Marker Heatmap
marker<-read.csv("Stage_marker.CSV")
head(marker)
targetcount<-normalized_counts[down,1:24]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
 
###sample corealation heatmap####
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)[,1:24]
#rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# ??????????pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
# ?㼶????
library(amap)
library(heatmap.2)
# ??ͼ????
pheatmap(pearson_cor,
         cluster_cols = F,cluster_rows = F,
         color = hmcol,
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)

#GSEA analysis
#Day2
MS_H1_Day2
gene=bitr(rownames(MS_H1_Day2),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=MS_H1_Day2$log2FoldChange, #??????foldchange
                      SYMBOL = rownames(MS_H1_Day2)) #??ס???Ļ?????ͷ????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID #ʹ??ת???õ?ID
geneList=sort(geneList,decreasing = T) #?Ӹߵ???????
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") #??gmt?ļ?
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.65,
           TERM2GENE = kegmt) #GSEA????
#dotplot(KEGG) #????ͼ 
#dotplot(KEGG,color="pvalue")  #??pֵ????ͼ 
#dotplot(KEGG,split=".sign")+facet_grid(~.sign) #????ͼ?????ҷ??漤????????
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #??????ʾ??ʽ
library(enrichplot)
#?ض?ͨ·??ͼ
pdf("Day2-GSEA-MS.pdf",width = 10,height = 10)
for (i in 1:13){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ????һ??????ά??ͼ??????ʾpֵ
dev.off()
write.csv(KEGG,"GSVA-KEGG-Day2.csv")

#Day5
MS_H1_Day5
gene=bitr(rownames(MS_H1_Day5),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=MS_H1_Day5$log2FoldChange, #??????foldchange
                      SYMBOL = rownames(MS_H1_Day5)) #??ס???Ļ?????ͷ????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID #ʹ??ת???õ?ID
geneList=sort(geneList,decreasing = T) #?Ӹߵ???????
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") #??gmt?ļ?
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.2,
           TERM2GENE = kegmt) #GSEA????
#dotplot(KEGG) #????ͼ 
#dotplot(KEGG,color="pvalue")  #??pֵ????ͼ 
#dotplot(KEGG,split=".sign")+facet_grid(~.sign) #????ͼ?????ҷ??漤????????
pdf("Day5-GSEA-MS-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #??????ʾ??ʽ
dev.off()
pdf("Day5-GSEA-MS.pdf",width = 10,height = 10)
for (i in 1:45){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ????һ??????ά??ͼ??????ʾpֵ
dev.off()
write.csv(KEGG,"GSVA-KEGG-Day5.csv")

#Day0
MS_H1_Day0
gene=bitr(rownames(MS_H1_Day0),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=MS_H1_Day0$log2FoldChange, #??????foldchange
                      SYMBOL = rownames(MS_H1_Day0)) #??ס???Ļ?????ͷ????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID #ʹ??ת???õ?ID
geneList=sort(geneList,decreasing = T) #?Ӹߵ???????
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") #??gmt?ļ?
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.9,
           TERM2GENE = kegmt) #GSEA????
write.csv(KEGG,"GSVA-KEGG-Day0.csv")

#dotplot(KEGG) #????ͼ 
#dotplot(KEGG,color="pvalue")  #??pֵ????ͼ 
#dotplot(KEGG,split=".sign")+facet_grid(~.sign) #????ͼ?????ҷ??漤????????
pdf("Day0-GSEA-MS-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #??????ʾ??ʽ
dev.off()
pdf("Day0-GSEA-MS.pdf",width = 10,height = 10)
for (i in 1:5){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ????һ??????ά??ͼ??????ʾpֵ
dev.off()

#Day10
MS_H1_Day10
gene=bitr(rownames(MS_H1_Day10),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=MS_H1_Day10$log2FoldChange, #??????foldchange
                      SYMBOL = rownames(MS_H1_Day10)) #??ס???Ļ?????ͷ????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID #ʹ??ת???õ?ID
geneList=sort(geneList,decreasing = T) #?Ӹߵ???????
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") #??gmt?ļ?
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.5,
           TERM2GENE = kegmt) #GSEA????
write.csv(KEGG,"GSVA-KEGG-Day10.csv")
#dotplot(KEGG) #????ͼ 
#dotplot(KEGG,color="pvalue")  #??pֵ????ͼ 
#dotplot(KEGG,split=".sign")+facet_grid(~.sign) #????ͼ?????ҷ??漤????????
pdf("Day10-GSEA-MS-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #??????ʾ??ʽ
dev.off()
pdf("Day10-GSEA-MS.pdf",width = 10,height = 10)
for (i in 1:53){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ????һ??????ά??ͼ??????ʾpֵ
dev.off()

#GSVA 
library(dplyr)
library(msigdbr)
library(GSVA)
library(patchwork)
m_df = msigdbr(species = "Homo sapiens", 
               category = "C5",
               subcategory = "GO:BP") #ѡȡ????????
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
expr=as.matrix(normalized_counts)
c5_BP <- gsva(expr, 
             msigdbr_list, 
             kcdf="Gaussian",
             method = "gsva",
             parallel.sz=10)
pheatmap(c5_BP[1768:1769,1:24],
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         color =colorRampPalette(c("blue", "white","red"))(100),
         cellwidth = 10, cellheight = 15,
         fontsize = 10)
grep("GOBP_HEART_DEVELOPMENT",rownames(c5_BP))



library(limma)
library(DESeq2)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
counts<-read.table("./TAD_all_rawcounts.txt")
head(counts)
# create pesudo replicates for batch2#
counts$H1_Day5_rep2_batch2<-counts$H1_Day5_rep1_batch2
counts$MS_Day5_rep2_batch2<-counts$MS_Day5_rep1_batch2
str(counts)
counts<-counts[,c(1:24,25,39,26,40,27:38)]
counts<-counts[,1:28]


###time course####
head(counts)
library(ImpulseDE2)

data = counts[apply(counts, 1, sum) > 1 , ] 
Sample = colnames(data)
Condition = c(rep("control",12),rep("case",12),rep("control",2),rep("case",2))
Time =c(rep(0,3),rep(2,3),rep(5,3),rep(10,3),
        rep(0,3),rep(2,3),rep(5,3),rep(10,3),rep(5,4))
Batch = c(rep("batch3",12),rep("batch1",12),rep("batch2",4))
dfAnnotation=data.frame(Sample,Condition,Time,Batch,row.names = Sample)
data = as.matrix(data)
###For transient 
objectImpulseDE5 <- runImpulseDE2(
  matCountData    = data, 
  dfAnnotation    = dfAnnotation,
  boolCaseCtrl    = TRUE,
  boolIdentifyTransients = TRUE,
  vecConfounders  = c("Batch"),
  scaQThres=0.05,
  scaNProc        = 12 )

p<-plotHeatmap(
  objectImpulseDE2       = objectImpulseDE5,
  strCondition           = "control",
  boolIdentifyTransients = T,
  scaQThres              = 0.01)

gene<-p$lsvecGeneGroups
normalized_counts<-read.table("normalized_counts_MS_H1.txt")
all<-c(gene$transition_up,gene$transition_down,gene$transient_up,gene$transient_down)
targetcount<-normalized_counts[all,c(1:24)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
condition <- factor(c(rep("H1",12),rep("MEIS1 KO",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))

annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)
rownames(annotation_col) = factor(colnames(count))
#rownames(df) = factor(df$gene)
type<-factor(c(rep("transition_up",2221),rep("transition_down",1999),
    rep("transient_up",1411),rep("transient_down",1054)))
annotation_row = data.frame(type=type)
rownames(annotation_row) = factor(all)
ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"),
  Group = c("H1"="#F39B7FB2","MEIS1 KO" ="#8491B4B2"),
  type=c("transition_up"="#FFAD60","transition_down"="#96CEB4",
         "transient_up"="#D9534F","transient_down"="#FFEEAD")
  )

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("Stage-heatmap-MEIS1-test.pdf",width = 6,height = 10)
pheatmap(count,
         cluster_cols = F,cluster_rows = F,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=F,show_colnames=T)
dev.off()




#3D PCA#
pca_count<-read.csv("../PCA/RNAseq-PCAcount.csv",row.names = 1)
library("FactoMineR")
library("factoextra")
library("scatterplot3d")
library("gmodels")
pca_count<-pca_count[,1:24]

pca.info <- fast.prcomp(pca_count)
head(pca.info$rotation)
head(pca.info)
pca.data <- data.frame(sample =rownames(pca.info$rotation),
                       Group = condition, 
                       Stage = timepoint,
                       pca.info$rotation)

table(pca.data$Group,pca.data$Stage)
library(ggsci)
library("scales")
colors=pal_npg("nrc")(10)
show_col(pal_npg("nrc")(10))
colors_pal<-colors[c(3,2,5,1)]
colors <- colors_pal[as.factor(pca.data$Stage)]

shape.lib=c(16,17,15)
shapes <- shape.lib[as.factor(pca.data$Group)]
 #计算PC值，并替换列名，用来替换坐标轴上的标签
 pVar <- pca.info$sdev^2/sum(pca.info$sdev^2)
 pVar = round(pVar,digits = 3)

   paste0("PC1 (",as.character(pVar[1] * 100 ),"%)")
   paste0("PC2 (",as.character(pVar[2] * 100 ),"%)")
   paste0("PC3 (",as.character(pVar[3] * 100 ),"%)")

str(pca.data)
s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],
                     pch =shapes, color = colors,mar=c(5,5,5,5),
                     angle = 60, type="p",cex.symbols = 1,
                     main = "MEIS1 PCA plot",
                     xlab="PC1 : 94.6% variance",
                     ylab = "PC2 : 2.7% variance",
                     zlab = "PC3 : 1.1% variance") 
legend("topright", legend = c("Day0","Day2","Day5","Day10"),
       col =colors_pal, pch = 15,
       inset = -0.2,xpd = TRUE, horiz = FALSE)
legend("bottomright", legend = c("H1","MEIS1 KO"),
       col ="black", pch = shape.lib,
       inset = -0.2,xpd = TRUE, horiz = FALSE)

lines3d(pca.data[1:4,c("PC1","PC2","PC3")])


points3d()

library(threejs)

## Not run: 
x <- rnorm(5)
y <- rnorm(5)
z <- rnorm(5)
scatterplot3js(x, y, z, pch="@", ) %>%
  lines3d(c(1, 2), c(3, 4), lwd=2)

scatterplot3js(pca.data[,"PC1"], 
               pca.data[,"PC2"], 
               pca.data[,"PC2"], pch =c("@","."), color = colors, ) %>%
  lines3d(c(1, 3), c(3, 4), lwd=2)

?lines3d
str(pca.data[,c("PC1","PC2","PC3")])
x
pca.data[,"PC1"]