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
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
dds <- DESeq(dds)
res <- results(dds)
setwd("D:/project/TAD/RNAseq/SA-Control")
write.csv(normalized_counts,"SA-control-normlized-counts.csv")
#condition DEG#
SA_Control_Day2<-results(dds, contrast=list( c("condition_SA_vs_Control",
                                               "timepointDay2.conditionSA") ))
SA_Control_Day2_diff_gene <-rownames(subset(SA_Control_Day2, 
                                            padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(SA_Control_Day2,"SA_Control_Day2.csv")
#volcano plot#
pdf("SA-vs-Control-Day2-volcano.pdf")
data<-as.data.frame(SA_Control_Day2)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'SA-KO','Control'),
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
       title=paste("Day2:SA-KO vs Control DEG","\n","Control:1110 genes; SA-KO:2013 genes;",sep=""))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%SA_Control_Day2_diff_gene),]
targetcount<-targetcount[,c(4:6,18:20)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("Day2-SA-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()

###Day0 
SA_Control_Day0<-results(dds, contrast=list( c("condition_SA_vs_Control") ))
SA_Control_Day0_diff_gene <-rownames(subset(SA_Control_Day0, 
                                            padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(SA_Control_Day0,"SA_Control_Day0.csv")

#volcano plot#
pdf("SA-vs-Control-Day0-volcano.pdf")
data<-as.data.frame(SA_Control_Day0)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'SA-KO','Control'),
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
       title=paste("Day0:SA-KO vs Control DEG","\n","Control:506 genes; SA-KO:1000 genes;",sep=""))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

targetcount<-normalized_counts[which(rownames(normalized_counts)%in%SA_Control_Day0_diff_gene),]
targetcount<-targetcount[,c(1:3,15:17)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("Day0-SA-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()

#Day5
SA_Control_Day5<-results(dds, contrast=list( c("condition_SA_vs_Control",
                                               "timepointDay5.conditionSA") ))
SA_Control_Day5_diff_gene <-rownames(subset(SA_Control_Day5, 
                                            padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(SA_Control_Day5,"SA_Control_Day5.csv")

#volcano plot#
pdf("SA-vs-Control-Day5-volcano.pdf")
data<-as.data.frame(SA_Control_Day5)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'SA-KO','Control'),
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
       title=paste("Day5:SA-KO vs Control DEG","\n","Control:1707 genes; SA-KO:1780 genes;",sep=""))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%SA_Control_Day5_diff_gene),]
targetcount<-targetcount[,c(7:9,21:23)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("Day5-SA-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()

#Day10:
SA_Control_Day10<-results(dds, contrast=list( c("condition_SA_vs_Control",
                                                "timepointDay10.conditionSA") ))
SA_Control_Day10_diff_gene <-rownames(subset(SA_Control_Day10, 
                                             padj < 0.05 & abs(log2FoldChange)> 1))
write.csv(SA_Control_Day10,"SA_Control_Day10.csv")

#volcano plot#
pdf("SA-vs-Control-Day10-volcano.pdf")
data<-as.data.frame(SA_Control_Day10)
head(data)
data<-na.omit(data)
data$change = ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'SA-KO','Control'),
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
       title=paste("Day10:SA-KO vs Control DEG","\n","Control:2165 genes; SA-KO:1900 genes;",sep=""))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Heatmap for DEG and GO KEGG #
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%SA_Control_Day10_diff_gene),]
targetcount<-targetcount[,c(10:12,24:26)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("Day10-SA-diff-logFC1-heatmap.pdf",width=6,height=10)
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
####GO and KEGG for SA-KO-DEG-Day2:
SA_Control_Day5_up <-subset(SA_Control_Day5, padj < 0.05 & log2FoldChange >1) 
SA_Control_Day5_up<-rownames(SA_Control_Day5_up)
SA_Control_Day5_up
pdf("SA_Control_Day5_up-GO.pdf")

gene.df <- bitr(SA_Control_Day5_up, fromType = "SYMBOL",
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
write.csv(ego,"SA_Control_Day5_upgene-BP.csv")


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
write.csv(ego,"SA_Control_Day5_upgene-MF.csv")

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
write.csv(ego,"SA_Control_Day5_upgene-CC.csv")
dev.off()
######KEGG##########
pdf("SA_Control_Day5_upgene-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"SA_Control_Day5_upgene-KEGG.csv")

dev.off()
########downgene############
SA_Control_Day5_down <-subset(SA_Control_Day5, padj < 0.05 & log2FoldChange < -1) 
str(SA_Control_Day5_down)
SA_Control_Day5_down<-rownames(SA_Control_Day5_down)
SA_Control_Day5_down:174 genes
pdf("SA_Control_Day5_downgene-GO.pdf")
gene.df <- bitr(SA_Control_Day5_down, fromType = "SYMBOL",
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
write.csv(ego,"SA_Control_Day5_downgene-GO-BP.csv")


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
write.csv(ego,"SA_Control_Day5-downgene-GO-MF.csv")

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
write.csv(ego,"SA_Control_Day5-downgene-GO-CC.csv")
dev.off()
######KEGG##########
pdf("SA_Control_Day5-downgene-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"SA_Control_Day5-downgene-KEGG.csv")
dev.off()


#time point DEG#
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
targetcount<-targetcount[,c(1:12,15:26)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
pdf("SA-diff-logFC1-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
p<-pheatmap(count,cluster_cols = F,cluster_rows = T,
            cutree_rows =6, 
            color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
            #cellwidth = 10, cellheight = 10,
            show_rownames=F,show_colnames=T)
row_cluster <- cutree(p$tree_row, k=6)
head(row_cluster)
Day0<-names(row_cluster[row_cluster==1])
Day2<-names(row_cluster[row_cluster==4])
Day5<-names(row_cluster[row_cluster==3])
Day10<-names(row_cluster[row_cluster==2])


count=t(scale(t(targetcount[c(Day0,Day2,Day5,Day10),]),scale = T,center = T))

pheatmap(count,cluster_cols = F,cluster_rows = F,
         cutree_rows =1, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
####GO and KEGG for SA-KO-Day0:
####GO and KEGG for SA-KO-Day0:
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
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"Day0-KEGG.csv")

dev.off()

write.csv(normalized_counts[,c(1:12,15:26)],"SA-Control-normalized_counts-cellNet.csv")


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
pdf("SA-Control-PCA.pdf",width = 8,height = 5)
p+theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()
head(assay(vsd))
pcacount<-assay(vsd)
head(pcacount)
write.table(pcacount,"SA-Control-RNAseq-PCAcount.csv")

#cell Net analysis
data<-read.csv("./CellNet Analysis/classification_results_2022-03-23.csv",row.names = 1)
str(data)
data<-data[,1:24]
condition <- factor(c(rep("Control",12),rep("SALL1 KD",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))

annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)
rownames(annotation_col) = factor(colnames(data))
ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"),
  Group = c("Control"="#F39B7FB2","SALL1 KD" ="#8491B4B2"))

pdf("SA-Control-CellNet.pdf",width = 10,height = 6)
pheatmap(data[c(1:12,14:15),],cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("black", "#48AF3A", "#E7E639"))(100),
         cellwidth = 15, cellheight = 10,
         annotation_col = annotation_col,
         annotation_colors = my47colors,
         #annotation_row = annotation_row,
         show_rownames=T,show_colnames=T)
dev.off()
ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"))


#Stage Marker Heatmap
marker<-read.csv("Stage_marker.CSV")
head(marker)
targetcount<-normalized_counts[marker$marker,c(1:12,15:26)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)
rownames(annotation_col) = factor(colnames(count))

annotation_row = data.frame(
  GeneClass = factor(marker$Stage)
)
rownames(annotation_row) = factor(marker$marker)
ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"),
  Group = c("Control"="#F39B7FB2","SALL1 KD" ="#8491B4B2"),
  GeneClass=c("Cardiac progenitor"="#00A087B2","cardiomyocytes" ="#E64B35B2",
              "Mesoderm"="#4DBBD5B2" ,"Pluripotency"="#3C5488B2"))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

pheatmap(count,cluster_cols = F,cluster_rows = F,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=T,show_colnames=T)

###SA-Neuron marker heatmap 
Ectoderm <- c("PAX6","OTX2","SOX1","FGF5","LHX2","OTX1","TUBB3","ENO2","MAP2","NEUROD1","NEUROG2","ZNF521",
              "MSI1","DCX","CRABP1","ALDH1A3","SOX3","OLIG1","CYP26A1","SIX3","FEZF1","FEZF2","LMX1A","TFAP2A",
              "SIX6","FOXG1","AMBN","DLX3","DLX6","PAX3","PAX7","ALDH1A2","SOX10")
targetcount<-normalized_counts[Ectoderm,c(1:12,15:26)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
condition <- factor(c(rep("Control",12),rep("SALL1 KD",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))

annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)
rownames(annotation_col) = factor(colnames(count))

annotation_row = data.frame(
  GeneClass = factor(rep("Ectoderm",length(Ectoderm)))
)
rownames(annotation_row) = factor(Ectoderm)
ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"),
  Group = c("Control"="#F39B7FB2","SALL1 KD" ="#8491B4B2"),
  GeneClass=c("Ectoderm"="#FFD24C"))

count<-na.omit(count)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_rownames=T,show_colnames=F)
###sample corealation heatmap####
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)[,c(1:12,15:26)]
#rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
library(amap)

condition <- factor(c(rep("WT",12),rep("MEIS1 KO",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))

annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)

rownames(annotation_col) = factor(colnames(pearson_cor))

ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"),
  Group = c("WT"="#F39B7FB2","MEIS1 KO" ="#8491B4B2"))


pheatmap(pearson_cor,
         cluster_cols = F,cluster_rows = F,
         color = hmcol,
         border_color = NA,
         annotation_col = annotation_col, 
         annotation_row = annotation_col,
         annotation_colors = ann_colors,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)

#GSEA analysis
#Day2
SA_Control_Day2
gene=bitr(rownames(SA_Control_Day2),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=SA_Control_Day2$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day2)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID #??1??????a????o??????ID
geneList=sort(geneList,decreasing = T) #?䨮?????????̨?????D??
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") #????gmt??????t
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.4,
           TERM2GENE = kegmt) 
write.csv(KEGG,"GSVA-KEGG-Day2.csv")
pdf("Day2-GSEA-SA-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)

pdf("Day2-GSEA-SA.pdf",width = 10,height = 10)
for (i in 1:21){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) 
dev.off()

#Day10
SA_Control_Day10
gene=bitr(rownames(SA_Control_Day10),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=SA_Control_Day10$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day10)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") 
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.4,
           TERM2GENE = kegmt) 
write.csv(KEGG,"GSVA-KEGG-Day10.csv")
pdf("Day10-GSEA-SA-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") 
dev.off()
library(enrichplot)

pdf("Day10-GSEA-SA.pdf",width = 10,height = 10)
for (i in 1:34){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) 
dev.off()
#Day5
SA_Control_Day5
gene=bitr(rownames(SA_Control_Day5),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=SA_Control_Day5$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day5)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") 
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.38,
           TERM2GENE = kegmt) 
write.csv(KEGG,"GSVA-KEGG-Day5.csv")
pdf("Day5-GSEA-SA-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") 
dev.off()
library(enrichplot)

pdf("Day5-GSEA-SA.pdf",width = 10,height = 10)
for (i in 1:20){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) 
dev.off()

#Day0
SA_Control_Day0
gene=bitr(rownames(SA_Control_Day0),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 

gene_df <- data.frame(logFC=SA_Control_Day0$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day0)) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") 
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.6,
           TERM2GENE = kegmt) 
write.csv(KEGG,"GSVA-KEGG-Day0.csv")
pdf("Day0-GSEA-SA-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") 
dev.off()
library(enrichplot)

pdf("Day0-GSEA-SA.pdf",width = 10,height = 10)
for (i in 1:13){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) 
dev.off()

###time course####
head(counts)
library(ImpulseDE2)

head(counts)
library(ImpulseDE2)
data = counts
Sample = colnames(data)
Condition = c(rep("case",12),rep("control",14))
Time =c(rep(0,3),rep(2,3),rep(5,3),rep(10,3),
        rep(0,5),rep(2,3),rep(5,3),rep(10,3))
#Batch = c(rep("batch4",14),rep("batch2",12))
dfAnnotation=data.frame(Sample,Condition,Time,row.names = Sample)
data = as.matrix(data)
###For transient 
objectImpulseDE5 <- runImpulseDE2(
  matCountData    = data, 
  dfAnnotation    = dfAnnotation,
  boolCaseCtrl    = TRUE,
  boolIdentifyTransients = TRUE,
  #vecConfounders  = c("Batch"),
  scaQThres=0.05,
  scaNProc        = 12 )

plotHeatmap(
  objectImpulseDE2       = objectImpulseDE5,
  strCondition           = "case",
  boolIdentifyTransients = T,
  scaQThres              = 0.01)

df<-read.csv("lsvecGeneGroups.csv")
head(df)
gene<-df$gene
targetcount<-normalized_counts[gene,c(1:12,15:26)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)
rownames(annotation_col) = factor(colnames(count))
#rownames(df) = factor(df$gene)
annotation_row = data.frame(type=df$type)
rownames(annotation_row) = factor(df$gene)
ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"),
  Group = c("Control"="#F39B7FB2","SALL1 KD" ="#8491B4B2"),
  type=c("transition_up"="#FFAD60","transition_down"="#96CEB4",
         "transient_up"="#D9534F","transient_down"="#FFEEAD")
  )

count<-na.omit(count)
library(pheatmap)
#df<-df[,-1:-2]
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("Stage-heatmap.pdf",width = 6,height = 10)
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
library(limma)
vsd <- vst(dds)
assay(vsd) <- limma::removeBatchEffect(assay(vsd),
                                       c(colData$batch))





pca_count<-read.csv("./SA-Control/SA-Control-RNAseq-PCAcount.csv",row.names = 1)
library("FactoMineR")
library("factoextra")
library("scatterplot3d")
library("gmodels")
  str(pca_count)
pca_count<-pca_count[,c(1:12,15:26)]

pca.info <- fast.prcomp(pca_count)
head(pca.info$rotation)
head(pca.info)
condition <- factor(c(rep("Control",12),rep("SALL1 KD",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))


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
s3d <- scatterplot3d(pca.data[,c("PC2","PC1","PC3")],
                     pch =shapes, color = colors,
                     mar=c(5,5,5,5),
                     angle = 60, type="p",cex.symbols = 1,
                     main = "SALL1 PCA plot",
                     xlab="PC2 : 2.4% variance",
                     ylab = "PC1 : 95.7% variance",
                     zlab = "PC3 : 0.7% variance") 

legend("topright", legend = c("Day0","Day2","Day5","Day10"),
       col =colors_pal, pch = 15,
       inset = -0.2,xpd = TRUE, horiz = FALSE)
legend("bottomright", legend = c("Control","SALL1 KD"),
       col ="black", pch = shape.lib,
       inset = -0.2,xpd = TRUE, horiz = FALSE)

