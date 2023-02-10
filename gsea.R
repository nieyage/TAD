##SALL1 GSEA on specified term 
library(clusterProfiler)
library(org.Hs.eg.db)
#Day2
SA_Control_Day2

gene_df <- data.frame(logFC=SA_Control_Day2$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day2)) 
#gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
#kegmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") 
reactomegmt<-read.gmt("c2.cp.reactome.v7.5.1.symbols.gmt")
all_gmt<-read.gmt("c5.all.v7.5.1.symbols.gmt")
reactome<-GSEA(geneList,pvalueCutoff = 0.37,
               TERM2GENE = reactomegmt)
write.csv(reactome,"GSVA-reactome-Day2.csv")

all<-GSEA(geneList,pvalueCutoff = 0.37,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day2.csv")
pdf("Day2-GSEA-SA-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day2-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:19){
  print(i)
  p<-gseaplot2(all,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,11:15,pvalue_table = T) 
dev.off()

#Day5 GSEA 
SA_Control_Day5
gene_df <- data.frame(logFC=SA_Control_Day5$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day5)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
all<-GSEA(geneList,pvalueCutoff = 0.05,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day5.csv")
pdf("Day5-GSEA-SA-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day5-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:68){
  print(i)
  p<-gseaplot2(all,i,base_size = 13,pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,6:10,pvalue_table = T) 
gseaplot2(all,11:15,pvalue_table = T) 
dev.off()
?gseaplot2
head(all_gmt)
GSEA(geneList,pvalueCutoff = 0.05,TERM2GENE = all_gmt)
#plot the specfied term GSEA 
gene_df <- data.frame(logFC=SA_Control_Day10$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day10)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
pdf("Day10-GSEA-SA-ectoderm-gseplot.pdf",width = 8,height = 6)
ectoderm<-all_gmt[grep("GOBP_ECTODERM_DEVELOPMENT",all_gmt$term),]
ectoderm_gsa<-GSEA(geneList,pvalueCutoff = 1,TERM2GENE = ectoderm)
gseaplot2(ectoderm_gsa,1,color="black",subplots =1:2,base_size = 10,pvalue_table = T)
gseaplot2(ectoderm_gsa,1,color="black",subplots =1:3,base_size = 10,pvalue_table = T)
dev.off()

unique(grep("EMBRYO",all_gmt$term,value = T))


#Day10
gene_df <- data.frame(logFC=SA_Control_Day10$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day10)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
all_gmt<-read.gmt("./c5.all.v7.5.1.symbols.gmt")
all<-GSEA(geneList,pvalueCutoff = 0.05,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day10.csv")
pdf("Day10-GSEA-SA-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day10-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:50){
  print(i)
  p<-gseaplot2(all,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,11:15,pvalue_table = T) 
dev.off()
#Day0
gene_df <- data.frame(logFC=SA_Control_Day0$log2FoldChange, 
                      SYMBOL = rownames(SA_Control_Day0)) 
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL 
geneList=sort(geneList,decreasing = T)
all_gmt<-read.gmt("../c5.all.v7.5.1.symbols.gmt")
all<-GSEA(geneList,pvalueCutoff = 0.5,TERM2GENE = all_gmt)
write.csv(all,"GSVA-all-Day0.csv")
pdf("Day0-GSEA-SA-dotplot-all.pdf",width = 16,height = 8)
dotplot(all,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
library(enrichplot)
pdf("Day0-GSEA-SA-all-gseplot.pdf",width = 10,height = 10)
for (i in 1:27){
  print(i)
  p<-gseaplot2(all,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(all,1:5,color="red",pvalue_table = T) 
gseaplot2(all,11:15,pvalue_table = T) 
dev.off()

