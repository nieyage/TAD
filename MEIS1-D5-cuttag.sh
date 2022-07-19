###extract the MEIS1 peak 

findMotifsGenome.pl MEIS1-D5_peaks.tmp hg19 ./ -find ../5_Motif/MEIS1-D5_motifDir_Homer/knownResults/known11.motif > MEIS1-D5_peaks_homer.txt;

awk '{print $1 }' MEIS1-D5_peaks_homer.txt |  sort -k 1 |uniq|grep -v "PositionID" | while read id;
do    read_id=$id"$";
    grep "$read_id" MEIS1-D5_peaks.tmp >> MEIS1-D5_peaks_homer.bed;done;



perl fragment_length_dist.pl HA-D5_sorted_rmDup_mapped_rmbl.bam MEIS1.fragL.txt
sort -n MEIS1.fragL.txt | uniq -c > MEIS1.frag.sort.txt

Rscript fragment_length_dist.R MEIS1.fragL.txt MEIS1.frag.sort.txt MEIS1.fragment_length_distribution.png MEIS1.fragment_length_distribution.txt
done

####
###MEIS1 Cut&tag 
awk '{print $1"\t"$2"\t"$3"\t"$4}' HA-D5_peaks.narrowPeak > MEIS1-D5_peaks.tmp

#bed to bw
sort -k1,1V -k2,2n -k3,3n MEIS1-D5_peaks_homer.bed>MEIS1-D5_peaks_homer_sort.bed

bedtools genomecov -i a.bed -g mm10.chrom.sizes.bed -bg > XX.bedgraph


computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/chip/SA-Chip/rawdata/3.align/hg19_RefSeq.bed  \
-S ../3_bigwig/HA-D5_CPM_normalized.bw \
--skipZeros  -o MEIS1-D5.gz  \
--outFileSortedRegions MEIS1-D5-bedtools-out.bed

plotHeatmap -m MEIS1-D5.gz  -out MEIS1-D5-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m MEIS1-D5.gz  -out MEIS1-D5-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 


library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, downstream=1000)

peakAnno <- annotatePeak("HA-D5_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("MEIS1-D5_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of MEIS1-binding loci\nrelative to TSS")

dev.off()

write.csv(as.data.frame(peakAnno),"MEIS1-cuttag-peak.annotation.csv",row.names = F)

####Motif Enrichment Analusis

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ../4_MACS2/HA-D5_peaks.narrowPeak > MEIS1-D5_peaks_homer.tmp
findMotifsGenome.pl MEIS1-D5_peaks_homer.tmp  hg19 MEIS1-D5_motifDir_Homer -len 8,10,12  

findMotifsGenome.pl MEIS1-D5_peaks_homer.bed  hg19 MEIS1-D5_motifDir_Homer_rmbg -len 8,10,12  

## BETA 

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' ../4_MACS2/HA-D5_peaks.narrowPeak > HA-D2_peaks.tmp

#cp /public/home/nieyg/project/TAD/cut_tag/20220525-SALL1-CUTTAG/4_MACS2/HA-D2_summits.bed MEIS1-D5_summits.bed
sed -i -e '/NA/d' MS_H1_Day5.txt

## DEG file first line need to add a "#"

 
conda activate beta_chip

awk '{print $1"\t"$2"\t"$3}' MEIS1-D5_peaks_homer.bed > MEIS1-D5_peaks_BETA.bed
BETA plus \
    -p MEIS1-D5_peaks_BETA.bed -e MS_H1_Day5.txt    -k O \
    --info 1,3,6\
    -g hg19 \
    --gs /public/home/nieyg/reference/genome/hg19/hg19.fa \
    --pn 4000 \
    --gname2 \
    -n MEIS1-Day5 \
    --df 0.05 \
    --da 1 \
    -o BETA_MEIS1_rmbg_pvalue


##MEME 
bedtools getfasta -fi /public/home/nieyg/reference/genome/hg19/hg19.fa \
-bed ../4_MACS2/HA-D5_peaks.narrowPeak  -fo MEIS1-D5.fa

# meme denovo motif 
meme  \
MEIS1-D5.fa \
-oc meme_denovo_max12 \
-dna \
-p 12 \
-maxw 12 \
-mod zoops \
-nmotifs 30 \
-revcomp

meme-chip \
-meme-p 12 \
-oc meme-chip-maxw12 \
-meme-maxw 12 \
-db /public/home/nieyg/database/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme \
MEIS1-D5.fa





data<-read.csv("MEIS1-Day5_downtarget.txt",sep="\t")
data<-data[data$rankproduct<0.05,]
gene<-data$GeneSymbol;
gene<-as.character(gene);

library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)


library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak("../MEIS1-D5_peaks_homer.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

bindinggene<-peakAnno@anno@ elementMetadata@listData$SYMBOL
DEG<-read.csv("MS_H1_Day5.csv");
DEG<-DEG[DEG$padj<0.05,];
DEG<-na.omit(DEG)

up<-DEG[DEG$log2FoldChange>1,]
down<-DEG[DEG$log2FoldChange < -1,]

up<-as.character(up$X)
down<-as.character(down$X)
up_overlap<-intersect(bindinggene,up)
down_overlap<-intersect(bindinggene,down)

pdf("MEIS1-Day5-downregulatedtarget-GO.pdf")
gene.df <- bitr(down_overlap, fromType = "SYMBOL",
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
write.csv(ego,"MEIS1-Day5-downregulatedtarget-BP.csv")


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
write.csv(ego,"MEIS1-Day5-downregulatedtarget-MF.csv")

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
write.csv(ego,"MEIS1-Day5-downregulatedtarget-CC.csv")
dev.off()
######KEGG##########
pdf("MEIS1-Day5-downregulatedtarget_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MEIS1-Day5-downregulatedtarget-KEGG.csv")

dev.off()

pdf("MEIS1-Day5-upregulatedtarget-GO.pdf")
gene.df <- bitr(up_overlap, fromType = "SYMBOL",
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
write.csv(ego,"MEIS1-Day5-upregulatedtarget-BP.csv")


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
write.csv(ego,"MEIS1-Day5-upregulatedtarget-MF.csv")

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
write.csv(ego,"MEIS1-Day5-upregulatedtarget-CC.csv")
dev.off()
######KEGG##########
pdf("MEIS1-Day5-upregulatedtarget_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MEIS1-Day5-upregulatedtarget-KEGG.csv")

dev.off()




vennplot<-venn.diagram(
  x = list(up,down,bindinggene),
  category.names = c("MEIS1 KO UP" , "MEIS1 KO DOWN" , "MEIS1 BINDING"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:3],
  alpha = 0.50,
  output=TRUE
)

pdf("vennplot.pdf")
grid.draw(vennplot)
dev.off()