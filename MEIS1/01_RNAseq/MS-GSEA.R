#MEIS1 GSEA 
#Day0 
library(clusterProfiler)
library(org.Hs.eg.db)
gene_df <- data.frame(logFC=MS_H1_Day0$log2FoldChange, 
                      SYMBOL = rownames(MS_H1_Day0)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
all_gmt<-read.gmt("../c5.all.v7.5.1.symbols.gmt")
all<-GSEA(geneList,pvalueCutoff = 0.5,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day0.csv")
pdf("Day0-GSEA-MS-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day0-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:19){
  print(i)
  p<-gseaplot2(all,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,11:15,pvalue_table = T) 
dev.off()

#Day2 
gene_df <- data.frame(logFC=MS_H1_Day2$log2FoldChange, 
                      SYMBOL = rownames(MS_H1_Day2)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
all<-GSEA(geneList,pvalueCutoff = 0.5,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day2.csv")
pdf("Day2-GSEA-MS-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day2-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:50){
  print(i)
  p<-gseaplot2(all,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,6:10,pvalue_table = T) 
dev.off()
gene_df <- data.frame(logFC=MS_H1_Day5$log2FoldChange, 
                      SYMBOL = rownames(MS_H1_Day5)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
all_gmt<-read.gmt("../c5.all.v7.5.1.symbols.gmt")
all<-GSEA(geneList,pvalueCutoff = 0.05,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day5.csv")
pdf("Day5-GSEA-MS-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day5-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:127){
  print(i)
  p<-gseaplot2(all,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,11:15,pvalue_table = T) 
dev.off()

#Day10
gene_df <- data.frame(logFC=MS_H1_Day10$log2FoldChange, 
                      SYMBOL = rownames(MS_H1_Day10)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
all_gmt<-read.gmt("../c5.all.v7.5.1.symbols.gmt")
all<-GSEA(geneList,pvalueCutoff = 0.05,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day10.csv")
pdf("Day10-GSEA-MS-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day10-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:30){
  print(i)
  p<-gseaplot2(all,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,11:15,pvalue_table = T) 
dev.off()

