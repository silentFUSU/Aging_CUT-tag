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
diff1 <- read.csv("data/samples/RNA/MEF_mid_age/diff_expression_gene_change_filter_bar.csv")
diff2 <- read.csv("data/samples/RNA/MEF_EZH2_inhibit/diff_expression_gene_change_filter_bar.csv")

diff1_increased <- diff1[which(diff1$Significant=="Up"),]
diff2_increased <- diff2[which(diff2$Significant=="Up"),]
common_increased <- intersect(diff1_increased$X,diff2_increased$X)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
common_increased <- bitr(common_increased,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
common_increased_GO <- enrichGO(common_increased$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)


diff1_decreased <- diff1[which(diff1$Significant=="Down"),]
diff2_decreased <- diff2[which(diff2$Significant=="Down"),]
common_decreased <- intersect(diff1_decreased$X,diff2_decreased$X)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
common_decreased <- bitr(common_decreased,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
common_decreased_GO <- enrichGO(common_decreased$ENTREZID,#GO富集分析
                                OrgDb = GO_database,
                                keyType = "ENTREZID",#设定读取的gene ID类型
                                ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                pvalueCutoff = 0.05,#设定p值阈值
                                qvalueCutoff = 0.05,#设定q值阈值
                                readable = T)


to_plot_increased <- common_increased_GO@result
to_plot_increased <- to_plot_increased[order(to_plot_increased$p.adjust),]
to_plot_increased <- to_plot_increased[c(1:3),]
to_plot_increased$p.adjust <- -log10(to_plot_increased$p.adjust)
to_plot_increased$label <- paste0(to_plot_increased$ID," ",to_plot_increased$Description)
to_plot_increased <- to_plot_increased[,c("label","p.adjust")]
p <- ggplot(to_plot_increased, aes(x = p.adjust, y = reorder(label, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "P Adjust", y = "Label") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
ggsave("result/figures/MEF_MEF_EZH2_inhibit_common_increased_gene_GO.pdf",p,width =15,height = 8)

to_plot_decreased <- common_decreased_GO@result
to_plot_decreased <- to_plot_decreased[order(to_plot_decreased$p.adjust),]
to_plot_decreased <- to_plot_decreased[c(1:3),]
to_plot_decreased$p.adjust <- -log10(to_plot_decreased$p.adjust)
to_plot_decreased$label <- paste0(to_plot_decreased$ID," ",to_plot_decreased$Description)
to_plot_decreased <- to_plot_decreased[,c("label","p.adjust")]
p <- ggplot(to_plot_decreased, aes(x = p.adjust, y = reorder(label, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "P Adjust", y = "Label") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
ggsave("result/figures/MEF_MEF_EZH2_inhibit_common_decreased_gene_GO.pdf",p,width =15,height = 8)


