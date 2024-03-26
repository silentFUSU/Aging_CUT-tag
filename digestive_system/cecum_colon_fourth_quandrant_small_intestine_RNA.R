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
library(gg.gap)
library(scales)
library(colorspace)  
library(UpSetR)
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
cecum <- read.delim("data/samples/cecum/H3K27me3_H3K9me3_intersect/10kb_all_significant_fourth_quadrant.bed",header = F)
colon <- read.delim("data/samples/colon/H3K27me3_H3K9me3_intersect/10kb_all_significant_fourth_quadrant.bed",header = F)

cecum_grange_obj <- GRanges(seqnames = cecum$V1,   
                      ranges = IRanges(start = cecum$V2, end = cecum$V3)) 
cecum_peakAnno <- annotatePeak(cecum_grange_obj, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
cecum_peakAnno<-unique(as.data.frame(cecum_peakAnno))

colon_grange_obj <- GRanges(seqnames = colon$V1,   
                            ranges = IRanges(start = colon$V2, end = colon$V3)) 
colon_peakAnno <- annotatePeak(colon_grange_obj, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
colon_peakAnno<-unique(as.data.frame(colon_peakAnno))

RNA <- read.csv("data/public_data/GSE132040/Small_Intestine.csv")
list <- list(colon=colon$V4,cecum=cecum$V4,RNA_increase=RNA$gene[which(RNA$Significant=="Up")],RNA_decrease=RNA$gene[which(RNA$Significant=="Down")])
upset(fromList(list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      nintersects = NA, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2 # 文字标签的大小
)
