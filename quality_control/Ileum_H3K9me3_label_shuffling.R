rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)
library(ChIPseeker)
library(clusterProfiler)
tissue <- "ileum"
antibody <- "H3K9me3"
df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_label_shuffling_diff.csv"))
unstable <- df[which(df$Significant!="Stable"),]
pheatmap::pheatmap(df[which(df$Significant=="Up"),c(7:10)],scale = "row",show_rownames = F)
pheatmap::pheatmap(df[which(df$Significant=="Down"),c(7:10)],scale = "row",show_rownames = F)
pheatmap::pheatmap(df[which(df$Significant=="Up"),c(7:10)],show_rownames = F)
pheatmap::pheatmap(df[which(df$Significant=="Down"),c(7:10)],show_rownames = F)

sample_correlation <- cor(df[which(df$Significant=="Up"),c(7:10)])
pheatmap::pheatmap(sample_correlation)
sample_correlation <- cor(df[which(df$Significant=="Down"),c(7:10)])
pheatmap::pheatmap(sample_correlation)

sample_correlation <- cor(df[,c(7:10)])
pheatmap::pheatmap(sample_correlation)

GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
out_up <- df[which(df$Significant=="Up"),]
up_peak <- GRanges(seqnames = out_up$Chr,   
                   ranges = IRanges(start = out_up$Start, end = out_up$End))
up_peak_anno <- annotatePeak(up_peak, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Mm.eg.db")
up_peak_anno <- unique(as.data.frame(up_peak_anno))
genelist_up <- bitr(up_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_up_GO,label_format = 50)


out_down <- df[which(df$Significant=="Down"),]
down_peak <- GRanges(seqnames = out_down$Chr,   
                   ranges = IRanges(start = out_down$Start, end = out_down$End))
down_peak_anno <- annotatePeak(down_peak, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Mm.eg.db")
down_peak_anno <- unique(as.data.frame(down_peak_anno))
genelist_down <- bitr(down_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_down_GO)
