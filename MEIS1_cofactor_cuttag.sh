macs2 callpeak -t P1_sorted_rmDup_mapped_rmbl.bam \
               -c IgG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n PBX1-IG 
               --outdir /public/home/nieyg/project/TAD/MEIS1/cuttag/PBX-TBX20-CUTTAG/4_MACS2 
macs2 callpeak -t P3_sorted_rmDup_mapped_rmbl.bam \
               -c IgG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n PBX3-IG 
               --outdir /public/home/nieyg/project/TAD/MEIS1/cuttag/PBX-TBX20-CUTTAG/4_MACS2 
macs2 callpeak -t T2_sorted_rmDup_mapped_rmbl.bam \
               -c IgG_sorted_rmDup_mapped_rmbl.bam  \
               -f BAMPE \
               -g hs \
               -n TBX20-IG 
               --outdir /public/home/nieyg/project/TAD/MEIS1/cuttag/PBX-TBX20-CUTTAG/4_MACS2 
cd /public/home/nieyg/project/TAD/MEIS1/cuttag/PBX-TBX20-CUTTAG/4_MACS2 
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ../2_bam/PBX1-IG_peaks.narrowPeak  > PBX1-IG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ../2_bam/PBX3-IG_peaks.narrowPeak  > PBX3-IG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ../2_bam/TBX20-IG_peaks.narrowPeak > TBX20-IG_peaks_homer.tmp
findMotifsGenome.pl PBX1-IG_peaks_homer.tmp   hg19 PBX1-IG_motifDir_Homer -len 6,8,10,12  
findMotifsGenome.pl PBX3-IG_peaks_homer.tmp   hg19 PBX3-IG_motifDir_Homer -len 6,8,10,12  
findMotifsGenome.pl TBX20-IG_peaks_homer.tmp  hg19 TBX20-IG_motifDir_Homer -len 6,8,10,12  

computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/pipeline/Chipseq/hg19_RefSeq.bed  \
-S IgG_CPM_normalized.bw P1_CPM_normalized.bw P3_CPM_normalized.bw T2_CPM_normalized.bw \
--skipZeros  -o MEIS1-cofactor.gz  \
--outFileSortedRegions MEIS1-cofactor-bedtools-out.bed
plotHeatmap -m MEIS1-cofactor.gz  -out MEIS1-cofactor-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m MEIS1-cofactor.gz  -out MEIS1-cofactor-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 


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
