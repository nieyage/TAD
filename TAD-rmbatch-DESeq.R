######remove batch#####
####DE remove batch effect (DEseq)####
library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
#####ºŸ…ËH1”ÎMS0œ‡Õ¨######

batch <- rep(c(rep("batch1",3),rep("batch2",3)),4)
#treat <- factor(rep(c(rep("WT",3),rep("MS",3)),4))
treat <- factor(rep(c(rep("WT",3),rep("MS",3)),4))
treat[4:6]<-"WT"
group<-factor(c(rep("Day0",6),rep("Day2",6),rep("Day5",6),rep("Day10",6)))
count<-read.table("D:/project/TAD/H1_MS_genesymbol_rawcount.txt",header=T)
colData <- data.frame(row.names=colnames(count),group=group,treat=treat,batch=batch)
colData
countData <- count[apply(count,1,sum) > 1 , ] 
head(countData)

dds<-DESeqDataSetFromMatrix(countData,colData, formula(~group+treat+batch)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
# log trans results######
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
#pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
#hc <- hclust(t(rlogMat) )
?hclust,method="pearson"

library(ggplot2)
vsd <- vst(dds)
head(vsd)
pca_data <- plotPCA(vsd, intgroup=c("group","treat","batch"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, shape =treat,color=group)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

####remove batch####
library(limma)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), c(colData$batch))
pca <- plotPCA(vsd, intgroup=c("batch","group","treat"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca, "percentVar"))
p<-ggplot(pca, aes(PC1, PC2, shape =treat,color=group)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

adjusted_counts <- sva::ComBat_seq(as.matrix(last), batch = batch, group = type)

head(assay(vsd))
pcacount<-assay(vsd)
head(pcacount)
write.table(pcacount,"RNAseq-PCAcount.csv")






####SVA################
library(sva)
batch <- rep(c(rep("1",3),rep("2",3)),4)
batch
group<-c(rep("Day0",6),rep("Day2",6),rep("Day5",6),rep("Day10",6))
mod = model.matrix(~as.factor(batch))  #group‰∏∫ÂàÜÁªÑ‰ø°ÊÅØ„ÄÇÊ≠§Ê≠•Êìç‰Ωú‰∏∫Âª∫Ê®°„Ä?
mod
modcombat = model.matrix(~1, data = count)
modcombat
exp2 = ComBat(dat=count, batch=batch, mod=modcombat, par.prior=TRUE, ref.batch="1")
head(exp2)
