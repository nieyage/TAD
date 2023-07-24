# signal of H3k27ac reduction region 
computeMatrix reference-point  --referencePoint center  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/5_analysis/H1_enrich.bed  \
-S /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/3_bigwig/WTac_CPM_merged.bw /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/3_bigwig/MKac_CPM_merged.bw \
--skipZeros  -o MEIS1-merged-H3k27ac-down-region.gz  \
--outFileSortedRegions MEIS1-merged-H3k27ac-down-region-out.bed
plotHeatmap -m MEIS1-merged-H3k27ac-down-region.gz  -out MEIS1-merged-H3k27ac-down-region-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m MEIS1-merged-H3k27ac-down-region.gz  -out MEIS1-merged-H3k27ac-down-region-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

# signal of H3k27ac increase region 
computeMatrix reference-point  --referencePoint center  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/5_analysis/MKO_enrich.bed  \
-S /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/3_bigwig/WTac_CPM_merged.bw /public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/3_bigwig/MKac_CPM_merged.bw \
--skipZeros  -o MEIS1-merged-H3k27ac-up-region.gz  \
--outFileSortedRegions MEIS1-merged-H3k27ac-up-region-out.bed
plotHeatmap -m MEIS1-merged-H3k27ac-up-region.gz  -out MEIS1-merged-H3k27ac-up-region-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m MEIS1-merged-H3k27ac-up-region.gz  -out MEIS1-merged-H3k27ac-up-region-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

# DEG in RNAseq Day5
MS_H1_Day5<- read.csv("/public/home/nieyg/project/TAD/MEIS1/rnaseq/MS_H1_Day5.csv")
MS_H1_Day5<- MS_H1_Day5[MS_H1_Day5$padj<0.05,]
MS_H1_Day5<- na.omit(MS_H1_Day5)
MS_H1_Day5_up<- unique(as.character(MS_H1_Day5[MS_H1_Day5$log2FoldChange>1,]$X))
MS_H1_Day5_down<- unique(as.character(MS_H1_Day5[MS_H1_Day5$log2FoldChange< -1,]$X))

# gene related to decreased region in H3k27ac 
# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
# Load data
samplefiles <- list.files("/public/home/nieyg/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/5_analysis", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("H1_enrich_MEIS1_target","H1", "MKO")
# Assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Get annotation data frame
H1_peakAnno <- annotatePeak(samplefiles[[2]],tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
MKO_peakAnno <- annotatePeak(samplefiles[[3]],tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
H1_annot <- as.data.frame(H1_peakAnno@anno)
MKO_annot <- as.data.frame(MKO_peakAnno@anno)
MKO_down_k27ac<- unique(H1_annot$SYMBOL)
MKO_up_k27ac<- unique(MKO_annot$SYMBOL)

intersect(MS_H1_Day5_down,MKO_down_k27ac)
intersect(MS_H1_Day5_up,MKO_up_k27ac)
intersect(MS_H1_Day5_up,MKO_down_k27ac)
intersect(MS_H1_Day5_down,MKO_up_k27ac)





