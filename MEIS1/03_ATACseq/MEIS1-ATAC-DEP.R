# MEIS1 KO and WT ATACseq DEP by DiffBind 

setwd("/Users/fraya/Documents/project/TAD/MEIS1/ATACseq/")
library(DiffBind)
library(Tidyverse)
dbObj <- dba(sampleSheet="diffbind_sample.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

pdf("MEIS1-ATAC-PCA-sample-correlation.pdf")
dba.plotPCA(dbObj,  attributes=DBA_TREATMENT, label=DBA_ID)
plot(dbObj)
dev.off()


# Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_TREATMENT,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)
#  summary of results
dba.show(dbObj, bContrasts=T)
pdf("MEIS1-ATAC-PCA-DEP.pdf")
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_TREATMENT, label=DBA_ID)
dev.off()
pdf("MEIS1-ATAC-DEP-MAplot.pdf")
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
pvals <- dba.plotBox(dbObj)
dev.off()

pdf("MEIS1-ATAC-DEP-plotVolcano.pdf")
dba.plotVolcano(dbObj)
dev.off()

comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="MKO_vs_WT_deseq2.txt", sep="\t", quote=F, col.names = NA)
deseq.bed <- out[ which(out$FDR < 0.05), c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="MKO_vs_WT_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)
library(dplyr)
WT_enrich <- out %>%  filter(FDR < 0.05 & Fold < 0) %>%  select(seqnames, start, end)
write.table(WT_enrich, file="WT_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
MKO_enrich <- out %>%  filter(FDR < 0.05 & Fold > 0) %>%  select(seqnames, start, end)
write.table(MKO_enrich, file="MKO_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

pdf("./MEIS1-ATAC-DEP-heatmap.pdf")
 corvals <- dba.plotHeatmap(dbObj)
 hmap <- colorRampPalette(c("white","blue"))(n = 13)
 readscores <- dba.plotHeatmap(dbObj, contrast=1, correlations=FALSE,scale="row", colScheme = hmap)
dev.off()
pdf("./MEIS1-ATAC-DEP-profiles.pdf")
profiles <- dba.plotProfile(dbObj)
dba.plotProfile(profiles)
profiles <- dba.plotProfile(tamoxifen,merge=DBA_REPLICATE)
dba.plotProfile(profiles)
dev.off()

library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb,upstream=1000, downstream=1000)
peakAnno <- annotatePeak("WT_enriched.bed",tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
MEIS1_KO_down<- unique(as.data.frame(peakAnno)$SYMBOL)
pdf("MEIS1-D5_closed.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of MEIS1-KO closed loci\nrelative to TSS")
dev.off()
write.csv(as.data.frame(peakAnno),"MEIS1-D5_closed_peak_annotation.csv",row.names = F)
peakAnno <- annotatePeak("MKO_enriched.bed",tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("MEIS1-D5_open.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of MEIS1-KO open loci\nrelative to TSS")
dev.off()
write.csv(as.data.frame(peakAnno),"MEIS1-D5_open_peak_annotation.csv",row.names = F)
MEIS1_KO_up<- unique(as.data.frame(peakAnno)$SYMBOL)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
gene.df <- bitr(MEIS1_KO_down, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
pdf("MEIS1_KO_closed-GO.pdf",width=8,height = 12)
    ego <- enrichGO(gene.df$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Hs.eg.db,
                    ont = "all",
                    pAdjustMethod  = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE);
    barplot(ego, x = "GeneRatio", #默认参数（x和color可以根据eG里面的内容更改）
            showCategory =10,
            split="ONTOLOGY") + 
      facet_grid(ONTOLOGY~., scale='free') 
dev.off()
  write.table(ego,"MEIS1_KO_closed-GO.csv")
pdf("MEIS1_KO_closed_KEGG.pdf",width=8,height = 5)
  ego <- enrichKEGG(
    gene = gene.df$ENTREZID,
    keyType = "kegg",
    organism  = "hsa",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
    p<-barplot(ego,x = "GeneRatio", color = "pvalue",  showCategory=10);
    print(p);
    edox<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
    write.table(edox,"MEIS1_KO_closed_KEGG.csv")
  }
  dev.off()


gene.df <- bitr(MEIS1_KO_up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
pdf("MEIS1_KO_open-GO.pdf",width=8,height = 12)
    ego <- enrichGO(gene.df$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Hs.eg.db,
                    ont = "all",
                    pAdjustMethod  = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE);
    barplot(ego, x = "GeneRatio", #默认参数（x和color可以根据eG里面的内容更改）
            showCategory =10,
            split="ONTOLOGY") + 
      facet_grid(ONTOLOGY~., scale='free') 
dev.off()
  write.table(ego,"MEIS1_KO_open-GO.csv")
pdf("MEIS1_KO_open_KEGG.pdf",width=8,height = 5)
  ego <- enrichKEGG(
    gene = gene.df$ENTREZID,
    keyType = "kegg",
    organism  = "hsa",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
    p<-barplot(ego,x = "GeneRatio", color = "pvalue",  showCategory=10);
    print(p);
    edox<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
    write.table(edox,"MEIS1_KO_open_KEGG.csv")
  }
  dev.off()


# motif enrichment 

findMotifsGenome.pl WT_enriched.bed hg19 MEIS1-KO-down_motifDir_Homer -len 6,8,10,12  
findMotifsGenome.pl MKO_enriched.bed hg19 MEIS1-KO-up_motifDir_Homer -len 6,8,10,12  

