### MEIS1 target peak accessibility
/public/home/nieyg/project/TAD/MEIS1/cuttag/20220608-MEIS1-CUTTAG/7_rmbg/MEIS1-D5_peaks_homer.bed

computeMatrix reference-point  --referencePoint center  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/MEIS1/cuttag/20220608-MEIS1-CUTTAG/7_rmbg/MEIS1-D5_peaks_homer.bed  \
-S IgG_CPM_normalized.bw P1_CPM_normalized.bw P3_CPM_normalized.bw T2_CPM_normalized.bw \
--skipZeros  -o MEIS1-center-cofactor.gz  \
--outFileSortedRegions MEIS1-center-cofactor-out.bed
plotHeatmap -m MEIS1-center-cofactor.gz  -out MEIS1-center-cofactor-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m MEIS1-center-cofactor.gz  -out MEIS1-center-cofactor-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

# MEIS1 TBX20 PBX target gene overlap and annotation 
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb,upstream=1000, downstream=1000)
peakAnno <- annotatePeak("/public/home/nieyg/project/TAD/MEIS1/cuttag/20220608-MEIS1-CUTTAG/7_rmbg/MEIS1-D5_peaks_homer.bed",tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
MEIS1_target<- unique(as.data.frame(peakAnno)$SYMBOL)
peakAnno <- annotatePeak("/public/home/nieyg/project/TAD/MEIS1/cuttag/PBX-TBX20-CUTTAG/2_bam/TBX20-IG_peaks.narrowPeak",tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
TBX20_target<- unique(as.data.frame(peakAnno)$SYMBOL)
pdf("TBX20_target.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of TBX20-binding closed loci\nrelative to TSS")
dev.off()

peakAnno <- annotatePeak("/public/home/nieyg/project/TAD/MEIS1/cuttag/PBX-TBX20-CUTTAG/2_bam/PBX1-IG_peaks.narrowPeak",tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
PBX1_target<- unique(as.data.frame(peakAnno)$SYMBOL)
pdf("PBX1_target.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of PBX1-binding closed loci\nrelative to TSS")
dev.off()
peakAnno <- annotatePeak("/public/home/nieyg/project/TAD/MEIS1/cuttag/PBX-TBX20-CUTTAG/2_bam/PBX3-IG_peaks.narrowPeak",tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
PBX3_target<- unique(as.data.frame(peakAnno)$SYMBOL)
pdf("PBX3_target.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of PBX1-binding closed loci\nrelative to TSS")
dev.off()

write.csv(MEIS1_target,"MEIS1_target-gene.csv")
write.csv(TBX20_target,"TBX20_target-gene.csv")
write.csv(PBX1_target,"PBX1_target-gene.csv")
write.csv(PBX3_target,"PBX3_target-gene.csv")


library(RColorBrewer)
library(VennDiagram)
vennplot<-venn.diagram(
  x = list(na.omit(MEIS1_target),na.omit(PBX1_target),na.omit(PBX3_target)),
  category.names = c("MEIS1" , "PBX1" ,"PBX3"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:3],
  alpha = 0.50,
  output=TRUE
)

pdf("3_TF_vennplot.pdf")
grid.draw(vennplot)
dev.off()

# all overlap 
inter <- get.venn.partitions(list(na.omit(MEIS1_target),na.omit(PBX1_target),na.omit(PBX3_target),na.omit(TBX20_target)))
gene<- inter$..values..[[1]]
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
write.csv(gene,"3TF_overlap_target_gene.csv")
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

pdf("MEIS1-cofactor-target-overlap-GO.pdf")
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
write.csv(ego,"MEIS1-cofactor-target-overlap-BP.csv")

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
write.csv(ego,"MEIS1-cofactor-target-overlap-MF.csv")

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
write.csv(ego,"MEIS1-cofactor-target-overlap-CC.csv")
dev.off()

######KEGG##########
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.1,
  use_internal_data =T)
pdf("MEIS1-cofactor-target-overlap_KEGG.pdf")
ego
barplot(ego, showCategory=20)
write.table(ego,"MEIS1-cofactor-target-overlap-KEGG.csv")
dev.off()
