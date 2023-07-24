# MEIS1 cut tag H3k27ac downstream analysis #
computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/reference/ref_gene/hg19_RefSeq.bed  \
-S WTac-1_CPM_normalized.bw WTac-2_CPM_normalized.bw MKac-1_CPM_normalized.bw MKac-2_CPM_normalized.bw \
--skipZeros  -o MEIS1-H3K27ac-CPM.gz  \
--outFileSortedRegions MEIS1-H3K27ac-bedtools-out.bed
plotHeatmap -m MEIS1-H3K27ac-CPM.gz  -out MEIS1-H3K27ac-CPM-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m MEIS1-H3K27ac-CPM.gz  -out MEIS1-H3K27ac-CPM-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

library(DiffBind)
sampleinfo<-read.csv("/public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/5_analysis/MEIS1_input_sampleinfo.csv")
MEIS1 <- dba(sampleSheet=sampleinfo)
MEIS1 <- dba.count(MEIS1, bUseSummarizeOverlaps=TRUE);
pdf("H3K27ac-PCA.pdf")
dba.plotPCA(MEIS1,  attributes=DBA_CONDITION, label=DBA_ID)
plot(MEIS1)
dev.off()
#MEIS1 <- dba.normalize(MEIS1)

#Establishing a model design and contrast
MEIS1 <- dba.contrast(MEIS1,categories=DBA_CONDITION,minMembers = 2);
dbObj <- dba.analyze(MEIS1, method=DBA_ALL_METHODS)
pdf("./MEIS1-H3k27ac-DiffBind-results.pdf")
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
dev.off()

res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=0.05,fold=1)

# Write to file

out <- as.data.frame(res_deseq)
write.table(out, file="./H1_vs_MKO_deseq2_FDR05FC1.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)
library(dplyr)
H1_enrich <- out %>% filter(FDR < 0.05 & Fold > 1) %>%   select(seqnames, start, end)
write.table(H1_enrich, file="H1_enrich.bed", sep="\t", quote=F, row.names=F, col.names=F)
MKO_enrich <- out %>%   filter(FDR < 0.05 & Fold < -1) %>%   select(seqnames, start, end)
write.table(MKO_enrich, file="MKO_enrich.bed", sep="\t", quote=F, row.names=F, col.names=F)

#get the sample normlizecount for heatmap 
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=0.05,fold=1,bCounts=T)
out <- as.data.frame(res_deseq)
data<-out[,12:15];
library(pheatmap);

filtermat=t(scale(t(data),scale = T,center = F))
pdf("./MEIS1-H3k27ac-DEseq-heatmap.pdf",width=5,heigh=8)
pheatmap(filtermat, cluster_rows=TRUE,
         show_rownames=F,border_color = NA,clustering_method="centroid",
         color = colorRampPalette(colors = c("white","blue"))(100),
         cluster_cols=T,cutree_rows =1,main = "H1 vs MKO",
         show_colnames = T,fontsize_row = 7.5)
dev.off()

# annotate the peak and gene functional analysis 

# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
# Load data
samplefiles <- list.files("/public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/5_analysis", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("H1", "MKO")

# Assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Get annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)
pdf("./Anno_peak_distribution.pdf",width=6,heigh=6)
plotAnnoPie(peakAnnoList[["H1"]]);
plotAnnoPie(peakAnnoList[["MKO"]]);
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList, title="Distribution of DEP relative to TSS")
dev.off()
# Get annotation data frame
H1_peakAnno <- annotatePeak(samplefiles[[1]],tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
MKO_peakAnno <- annotatePeak(samplefiles[[2]],tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
H1_annot <- as.data.frame(H1_peakAnno@anno)
MKO_annot <- as.data.frame(MKO_peakAnno@anno)
write.csv(H1_annot,"H1-peak-anno.csv")
write.csv(MKO_annot,"MKO-peak-anno.csv")


pdf("H1-H3K27ac-GO.pdf")

ego <- enrichGO(H1_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"H1-H3K27ac-BP.csv")

ego <- enrichGO(H1_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"H1-H3K27ac-MF.csv")

ego <- enrichGO(H1_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"H1-H3K27ac-CC.csv")
dev.off()
######KEGG##########
pdf("H1-H3K27ac_KEGG.pdf")
ego <- enrichKEGG(
  gene = H1_annot$geneId,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"H1-H3K27ac-KEGG.csv")
dev.off()

pdf("MKO-H3K27ac-GO.pdf")
ego <- enrichGO(MKO_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"MKO-H3K27ac-BP.csv")

ego <- enrichGO(MKO_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"MKO-H3K27ac-MF.csv")

ego <- enrichGO(MKO_annot$geneId,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"MKO-H3K27ac-CC.csv")
dev.off()
######KEGG##########
pdf("MKO-H3K27ac_KEGG.pdf")

ego <- enrichKEGG(
  gene = MKO_annot$geneId,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MKO-H3K27ac-KEGG.csv")

dev.off();

## Create a list with genes from each sample
#genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#
## Run KEGG analysis
#compKEGG <- compareCluster(geneCluster = genes, 
#                         fun = "enrichKEGG",
#                         organism = "human",
#                         pvalueCutoff  = 0.05, 
#                         pAdjustMethod = "BH")
#dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

# Motif Enrichment 

findMotifsGenome.pl MKO_enrich.bed  hg19 MEIS1-KO_motifDir_Homer -len 8,10,12  
findMotifsGenome.pl H1_enrich.bed  hg19 H1_motifDir_Homer -len 8,10,12  

bedtools intersect -a H1_enrich.bed -b /public/home/nieyg/project/TAD/MEIS1/cuttag/20220608-MEIS1-CUTTAG/4_MACS2/MEIS1-D5_peaks_homer.bed > H1_enrich_MEIS1_target.bed

peakAnno <- annotatePeak("H1_enrich_MEIS1_target.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

H1_enrich_MEIS1_target_gene<-peakAnno@anno@ elementMetadata@listData$SYMBOL

# H1_enrich_MEIS1_target_gene RNAseq heatmap and GO KEGG 

MEIS1_target_gene <- annotatePeak("/public/home/nieyg/project/TAD/MEIS1/cuttag/20220608-MEIS1-CUTTAG/4_MACS2/MEIS1-D5_peaks_homer.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
MEIS1_target_gene <- MEIS1_target_gene@anno$SYMBOL
H1gene <- annotatePeak("H1_enrich.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
H1gene <- H1gene@anno$SYMBOL;

overlap<-intersect(H1gene,MEIS1_target_gene)
pdf("H1_enrich_MEIS1_target_gene-GO.pdf")
gene.df <- bitr(overlap, fromType = "SYMBOL",
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
write.csv(ego,"H1_enrich_MEIS1_target_gene-BP.csv")


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
write.csv(ego,"H1_enrich_MEIS1_target_gene-MF.csv")

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
write.csv(ego,"H1_enrich_MEIS1_target_gene-CC.csv")
dev.off()
######KEGG##########
pdf("H1_enrich_MEIS1_target_gene_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H1_enrich_MEIS1_target_gene-KEGG.csv")

dev.off()

# heatmap for overlap gene 
# overlap with MKO down gene 
Day5_DEG<- read.csv("/public/home/nieyg/project/TAD/MEIS1/rnaseq/MS_H1_Day5.csv",row.names=1,header=T)
MEIS1_KO_downgene<-rownames(Day5_DEG[which(Day5_DEG$log2FoldChange< -1&Day5_DEG$padj<0.05),])

overlap<-intersect(H1gene,MEIS1_target_gene)
last_overlap<-intersect(MEIS1_KO_downgene,overlap)
library(VennDiagram)
library(RColorBrewer)
vennplot<-venn.diagram(
  x = list(na.omit(H1gene),na.omit(MEIS1_target_gene),MEIS1_KO_downgene),
  category.names = c("MEIS1 KO H3k27ac down" , "MEIS1 target gene" , "MEIS1 KO downregulated genes"),
  filename =NULL,
  fill = brewer.pal(7, "Set3")[1:3],
  alpha = 0.50,
  output=TRUE
)

pdf("MEIS1-vennplot.pdf")
grid.draw(vennplot)
dev.off()



library(pheatmap)


normalized_counts<-read.table("/public/home/nieyg/project/TAD/MEIS1/rnaseq/normalized_counts_MS_H1.txt",header=T)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%last_overlap),]
targetcount<-targetcount[,c(7:9,19:21)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
condition <- factor(c(rep("H1",3),rep("MEIS1 KO",3)))

annotation_col = data.frame(
  Group = condition
)
rownames(annotation_col) = factor(colnames(count))
ann_colors= list(Group = c("H1"="#F39B7FB2","MEIS1 KO"="#8491B4B2"))

count<-na.omit(count)
library(pheatmap)
#bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf("cuttag-overlap-gene-heatmap.pdf",width=6,height=10)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         #legend_breaks=seq(-2,2,1),
         #breaks=bk,
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         #cutree_rows=2,
         annotation_colors = ann_colors,
         show_rownames=T,show_colnames=F)

dev.off()

#the last overlap gene GO and KEGG 
pdf("MEIS1_regulate_gene-GO.pdf")
gene.df <- bitr(last_overlap, fromType = "SYMBOL",
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
write.csv(ego,"MEIS1_regulate_gene-BP.csv")


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
write.csv(ego,"MEIS1_regulate_gene-MF.csv")

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
write.csv(ego,"MEIS1_regulate_gene-CC.csv")
dev.off()
######KEGG##########
pdf("MEIS1_regulate_gene_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.1)
ego
barplot(ego, showCategory=20)
write.csv(ego,"MEIS1_regulate_gene-KEGG.csv")
dev.off();

#sample gene track 
#cardiac muscle tissue development: BMP2/JPH2/SORBS2/FZD7/TNNI1/S1PR1/NKX2-6
#compare with input:
nohup macs2 callpeak -t /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/MKac-1_sorted_rmDup_mapped_rmbl.bam \
-c /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam -f BAMPE -g hs \
--keep-dup all -n MKac_rmbg-1 --outdir /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/4_MACS2/ &

nohup macs2 callpeak -t /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/MKac-2_sorted_rmDup_mapped_rmbl.bam \
-c /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam -f BAMPE -g hs \
--keep-dup all -n MKac_rmbg-2 --outdir /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/4_MACS2/ &

nohup macs2 callpeak -t /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/WTac-1_sorted_rmDup_mapped_rmbl.bam \
-c /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam -f BAMPE -g hs \
--keep-dup all -n WTac_rmbg-1 --outdir /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/4_MACS2/ &

nohup macs2 callpeak -t /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/WTac-2_sorted_rmDup_mapped_rmbl.bam \
-c /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam -f BAMPE -g hs \
--keep-dup all -n WTac_rmbg-2 --outdir /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/4_MACS2/ &

#idr for rep

idr --samples MKac_rmbg-1_peaks.narrowPeak MKac_rmbg-2_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file MKac_rmbg-idr \
--plot \
--log-output-file MKac_rmbg.idr.log 

idr --samples WTac_rmbg-1_peaks.narrowPeak WTac_rmbg-2_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file WTac_rmbg-idr \
--plot \
--log-output-file WTac_rmbg.idr.log

# merge replicated peak 
bedtools intersect \
-a MKac_rmbg-1_peaks.narrowPeak \
-b MKac_rmbg-2_peaks.narrowPeak \
-u > MKac_rmbg-overlaps.bed

bedtools intersect \
-a WTac_rmbg-1_peaks.narrowPeak \
-b WTac_rmbg-2_peaks.narrowPeak \
-u > WTac_rmbg-overlaps.bed

#merge bw 
bigWigMerge WTac-*_CPM_normalized.bw WTac_CPM_normalized.bedGraph
bigWigMerge MKac-*_CPM_normalized.bw MKac_CPM_normalized.bedGraph

sort -k1,1 -k2,2n  WTac_CPM_normalized.bedGraph >  WTac_CPM_normalized_sorted.bedgraph
sort -k1,1 -k2,2n  MKac_CPM_normalized.bedGraph >  MKac_CPM_normalized_sorted.bedgraph
bedGraphToBigWig WTac_CPM_normalized_sorted.bedgraph /public/home/nieyg/reference/chrom_sizes/hg19.chrom.sizes WTac_CPM_merged.bw
bedGraphToBigWig MKac_CPM_normalized_sorted.bedgraph /public/home/nieyg/reference/chrom_sizes/hg19.chrom.sizes MKac_CPM_merged.bw

computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/reference/ref_gene/hg19_RefSeq.bed  \
-S WTac_CPM_merged.bw MKac_CPM_merged.bw \
--skipZeros  -o MEIS1-merged.gz  \
--outFileSortedRegions MEIS1-H3K27ac-bedtools-out.bed
plotHeatmap -m MEIS1-merged.gz  -out MEIS1-merged-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m MEIS1-merged.gz  -out MEIS1-merged-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 


