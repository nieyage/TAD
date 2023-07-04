#BETA pipeline

# install BETA 
conda create -y -n beta_chip python=2.7.15
conda install -y -c hcc beta 
conda install -y libiconv

# 对应BETA的三种功能，软件主要有三种子命令：
# 
#     BETA Basic: BETA Basic和BETA plus需要2个输入文件:
#         输入文件一，TF结合位点，即peak的位置；
#         输入文件二，敲降或敲除或过表达或激活TF后所有基因的变化倍数和P值;
#         BETA Basic能预测转录因子的激活/抑制功能，并识别直接靶基因。
#     BETA plus: 除了Basic的两个功能外，还能鉴定转录因子motif及其collaborator，相当于Basic的扩展。
#     BETA minus: 只有ChIPseq的数据的情况下，BETA minus可以根据bed文件计算出TF对靶基因的调控潜力。

#For SALL1
input file:
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' /public/home/nieyg/project/TAD/cut_tag/20220525-SALL1-CUTTAG/4_MACS2/HA-D2_peaks.narrowPeak > HA-D2_peaks.tmp

#cp /public/home/nieyg/project/TAD/cut_tag/20220525-SALL1-CUTTAG/4_MACS2/HA-D2_summits.bed SALL1-D2_summits.bed
sed -i -e '/NA/d' SA_Control_Day2.txt

## DEG file first line need to add a "#"

 
conda activate beta_chip


BETA plus \
    -p HA-D2_peaks.tmp \
    -e SA_Control_Day2.txt  \
    -k O \
    --info 1,3,7\
    -g hg19 \
    --gs /public/home/nieyg/reference/genome/hg19/hg19.fa \
    --pn 77870 \
    --gname2 \
    -n SALL1-Day2 \
    --df 0.05 \
    -o BETA_SALL1_Day2 


BETA plus \
    -p ../HA-D2_peaks.tmp \
    -e ../SA_Control_Day2.txt  \
    -k O \
    --info 1,3,7\
    -g hg19 \
    --gs /public/home/nieyg/reference/genome/hg19/hg19.fa \
    --pn 77870 \
    --gname2 \
    -n SALL1-Day2-nosummit \
    -o BETA_output_dir_nosummit \

##MEME 
bedtools getfasta -fi /public/home/nieyg/reference/genome/hg19/hg19.fa \
-bed /public/home/nieyg/project/TAD/cut_tag/20220525-SALL1-CUTTAG/4_MACS2/HA-D2_peaks.narrowPeak \
 -fo SALL1-D2.fa


# meme denovo motif 
meme  \
SALL1-D2.fa \
-oc meme_denovo_max12 \
-dna \
-p 12 \
-maxw 12 \
-mod zoops \
-nmotifs 30 \
-revcomp

[-w <w>]		motif width
	[-minw <minw>]		minimum motif width
	[-maxw <maxw>]		maximum motif width


meme-chip \
-meme-p 12 \
-oc meme-chip-maxw12 \
-meme-maxw 12 \
-db /public/home/nieyg/database/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme \
SALL1-D2.fa

#input last bam 

perl fragment_length_dist.pl HA-D2_sorted_rmDup_mapped_rmbl.bam SALL1.fragL.txt
sort -n SALL1.fragL.txt | uniq -c > SALL1.frag.sort.txt

Rscript fragment_length_dist.R SALL1.fragL.txt SALL1.frag.sort.txt SALL1.fragment_length_distribution.png SALL1.fragment_length_distribution.txt
done


library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

data<-read.csv("SALL1-Day2_downtarget.txt",sep="\t")
data<-data[data$rankproduct<0.05,]
gene<-data$GeneSymbol[1:1000];
gene<-as.character(gene);



pdf("SALL1-Day2-downregulatedtarget-GO.pdf")
gene.df <- bitr(gene, fromType = "SYMBOL",
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
write.csv(ego,"SALL1-Day2-downregulatedtarget-BP.csv")

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
write.csv(ego,"SALL1-Day2-downregulatedtarget-MF.csv")

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
write.csv(ego,"SALL1-Day2-downregulatedtarget-CC.csv")
dev.off()
######KEGG##########
pdf("SALL1-Day2-downregulatedtarget_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"SALL1-Day2-downregulatedtarget-KEGG.csv")

dev.off()

data<-read.csv("SALL1-Day2_uptarget.txt",sep="\t")
data<-data[data$rankproduct<0.05,]
gene<-data$GeneSymbol[1:2000];
gene<-as.character(gene);
pdf("SALL1-Day2-upregulatedtarget-GO.pdf")
gene.df <- bitr(gene, fromType = "SYMBOL",
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
write.csv(ego,"SALL1-Day2-upregulatedtarget-BP.csv")

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
write.csv(ego,"SALL1-Day2-upregulatedtarget-MF.csv")

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
write.csv(ego,"SALL1-Day2-upregulatedtarget-CC.csv")
dev.off()
######KEGG##########
pdf("SALL1-Day2-upregulatedtarget_KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"SALL1-Day2-upregulatedtarget-KEGG.csv")

dev.off()

