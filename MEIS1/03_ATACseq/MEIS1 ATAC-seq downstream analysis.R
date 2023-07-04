# MEIS1 ATAC-seq downstream analysis #
library(DiffBind)
sampleinfo<-read.csv("/public/home/nieyg/project/TAD/MEIS1/atacseq/7_analysis/DiffBind/MEIS1_input_sampleinfo.csv")
MEIS1 <- dba(sampleSheet=sampleinfo)
MEIS1 <- dba.count(MEIS1, bUseSummarizeOverlaps=TRUE)
pdf("MEIS1-ATAC-PCA.pdf")
dba.plotPCA(MEIS1,  attributes=DBA_CONDITION, label=DBA_ID)
plot(MEIS1)
dev.off()


#Establishing a model design and contrast
MEIS1 <- dba.contrast(MEIS1,categories=DBA_CONDITION,minMembers = 2);
dbObj <- dba.analyze(MEIS1, method=DBA_ALL_METHODS)
pdf("MEIS1-ATAC-DiffBind-results.pdf")
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
dev.off()


