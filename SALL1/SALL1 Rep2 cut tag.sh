#SALL1 Rep2 cut tag 
#Step1: Basical Figure

computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/chip/SA-Chip/rawdata/3.align/hg19_RefSeq.bed  \
-S SA_CPM_normalized.bw \
--skipZeros  -o SA-CPM.gz  \
--outFileSortedRegions SA-CPM-bedtools-out.bed
plotHeatmap -m SA-CPM.gz  -out SA-CPM-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m SA-CPM.gz  -out SA-CPM-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/chip/SA-Chip/rawdata/3.align/hg19_RefSeq.bed  \
-S D2-IG_CPM_normalized.bw \
--skipZeros  -o D2-IG-CPM.gz  \
--outFileSortedRegions D2-IG-CPM-bedtools-out.bed
plotHeatmap -m D2-IG-CPM.gz  -out D2-IG-CPM-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m D2-IG-CPM.gz  -out D2-IG-CPM-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/chip/SA-Chip/rawdata/3.align/hg19_RefSeq.bed  \
-S SA_CPM_normalized.bw D2-IG_CPM_normalized.bw \
--skipZeros  -o SA-D2-IG-CPM.gz  \
--outFileSortedRegions SA-D2-IG-bedtools-out.bed
plotHeatmap -m SA-D2-IG-CPM.gz  -out SA-D2-IG-CPM-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m SA-D2-IG-CPM.gz  -out SA-D2-IG-CPM-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 


library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, downstream=1000)

peakAnno <- annotatePeak("SA-D2-IG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("SA-D2-IG_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of SALL1-binding loci\nrelative to TSS")

dev.off()

write.csv(as.data.frame(peakAnno),"SA-D2-IG-cuttag-peak.annotation.csv",row.names = F)
