rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(stringr)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
diff <- read.csv("data/samples/RNA/MEF/diff_expression_gene_change_filter_bar.csv")

diff_Up <- diff[which(diff$Significant=="Up"),]
diff_Down <- diff[which(diff$Significant=="Down"),]
diff_Up <- diff_Up[order(diff_Up$fdr),]
diff_Down <- diff_Down[order(diff_Down$fdr),]




genelist_up <- bitr(diff_Up$X[1:1000],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,label_format = 50,showCategory = 20)

genelist_down <- bitr(diff$X[which(diff$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
barplot(genelist_down_GO,label_format = 50,showCategory = 20)



to_plot_increased <- genelist_up_GO@result
to_plot_increased <- to_plot_increased[order(to_plot_increased$p.adjust,to_plot_increased$pvalue),]
to_plot_increased <- to_plot_increased[c(1:5),]
to_plot_increased$p.adjust <- -log10(to_plot_increased$p.adjust)
to_plot_increased$label <- paste0(to_plot_increased$ID," ",to_plot_increased$Description)
to_plot_increased <- to_plot_increased[,c("label","p.adjust")]
p <- ggplot(to_plot_increased, aes(x = p.adjust, y = reorder(label, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "P Adjust", y = "Label") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))

to_plot_decreased <- genelist_down_GO@result
to_plot_decreased <- to_plot_decreased[order(to_plot_decreased$p.adjust),]
to_plot_decreased <- to_plot_decreased[c(1:5),]
to_plot_decreased$p.adjust <- -log10(to_plot_decreased$p.adjust)
to_plot_decreased$label <- paste0(to_plot_decreased$ID," ",to_plot_decreased$Description)
to_plot_decreased <- to_plot_decreased[,c("label","p.adjust")]
p <- ggplot(to_plot_decreased, aes(x = p.adjust, y = reorder(label, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "P Adjust", y = "Label") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
