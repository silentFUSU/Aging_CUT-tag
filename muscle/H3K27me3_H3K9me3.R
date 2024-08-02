rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(UpSetR)
tissue <- "thymus"

H3K27me3 <-  read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
H3K27me3_increase_bar <- unique(H3K27me3$Geneid[which(H3K27me3$Significant_bar=="Up")])
H3K27me3_decrease_bar <- unique(H3K27me3$Geneid[which(H3K27me3$Significant_bar=="Down")])
H3K9me3_increase_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Up")])
H3K9me3_decrease_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Down")])
list <- list(H3K27me3_increase=H3K27me3_increase_bar,H3K27me3_decrease=H3K27me3_decrease_bar,
             H3K9me3_increase=H3K9me3_increase_bar,H3K9me3_decrease=H3K9me3_decrease_bar)
data_type <- c("H3K27me3_increase","H3K27me3_decrease","H3K9me3_increase","H3K9me3_decrease")
upset(fromList(list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      sets=data_type,nsets = 100,     # 绘制的最大集合个数
      nintersects = 20, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 1.8 # 文字标签的大小
)
H3K27me3_H3K9me3_down <- intersect(H3K27me3_decrease_bar,H3K9me3_decrease_bar)
H3K27me3_H3K9me3_up <- intersect(H3K27me3_increase_bar,H3K9me3_increase_bar)
H3K27me3_up_H3K9me3_down <- intersect(H3K27me3_increase_bar,H3K9me3_decrease_bar)
mm10_10kb <- read.table("~/ref_data/mm10_10kb_bins.bed")
H3K27me3_H3K9me3_down <- mm10_10kb[which(mm10_10kb$V4 %in% H3K27me3_H3K9me3_down),]
H3K27me3_H3K9me3_up <- mm10_10kb[which(mm10_10kb$V4 %in% H3K27me3_H3K9me3_up),]
H3K27me3_up_H3K9me3_down  <- mm10_10kb[which(mm10_10kb$V4 %in% H3K27me3_up_H3K9me3_down),]
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
H3K27me3_H3K9me3_down_grange <- GRanges(seqnames = H3K27me3_H3K9me3_down$V1,   
                                         ranges = IRanges(start = H3K27me3_H3K9me3_down$V2, end = H3K27me3_H3K9me3_down$V3))
H3K27me3_H3K9me3_down_peak_anno <- annotatePeak(H3K27me3_H3K9me3_down_grange, tssRegion=c(-3000, 3000),
                                                 TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoBar(H3K27me3_H3K9me3_down_peak_anno)

H3K27me3_H3K9me3_up_grange <- GRanges(seqnames = H3K27me3_H3K9me3_up$V1,   
                                       ranges = IRanges(start = H3K27me3_H3K9me3_up$V2, end = H3K27me3_H3K9me3_up$V3))
H3K27me3_H3K9me3_up_peak_anno <- annotatePeak(H3K27me3_H3K9me3_up_grange, tssRegion=c(-3000, 3000),
                                               TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoBar(H3K27me3_H3K9me3_up_peak_anno)

plotAnnoBar(list(H3K27me3_H3K9me3_up=H3K27me3_H3K9me3_up_peak_anno,
                 H3K27me3_H3K9me3_down=H3K27me3_H3K9me3_down_peak_anno))

H3K27me3_H3K9me3_down_peak_anno <- unique(as.data.frame(H3K27me3_H3K9me3_down_peak_anno))
genelist_down <- bitr(H3K27me3_H3K9me3_down_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                             OrgDb = GO_database,
                             keyType = "ENTREZID",#设定读取的gene ID类型
                             ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                             pvalueCutoff = 0.05,#设定p值阈值
                             qvalueCutoff = 0.05,#设定q值阈值
                             readable = T)
barplot(genelist_down_GO)
