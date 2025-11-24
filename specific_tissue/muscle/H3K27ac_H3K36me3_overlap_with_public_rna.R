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
H3K27ac <- read.csv("data/samples/muscle/H3K27ac/H3K27ac_macs_narrowpeak_diff_FCbar_after_remove_batch_effect.csv")
H3K27ac_increase_bar <- unique(H3K27ac$symbol[which(H3K27ac$Significant_bar=="Up")])
H3K27ac_decrease_bar <- unique(H3K27ac$symbol[which(H3K27ac$Significant_bar=="Down")])
RNA_increase <- out$gene[which(out$Significant=="Up")]
RNA_decrease <- out$gene[which(out$Significant=="Down")]

H3K36me3 <- read.csv("data/samples/muscle/H3K36me3/H3K36me3_merge-W1000-G3000-E100_diff_FCbar_after_remove_batch_effect.csv")
H3K36me3_increase_bar <- unique(H3K36me3$symbol[which(H3K36me3$Significant_bar=="Up")])
H3K36me3_decrease_bar <- unique(H3K36me3$symbol[which(H3K36me3$Significant_bar=="Down")])
H3K36me3_increase <- unique(H3K36me3$symbol[which(H3K36me3$Significant=="Up")])
H3K36me3_decrease <- unique(H3K36me3$symbol[which(H3K36me3$Significant=="Down")])

increase_upset <- list(H3K36me3_increase,RNA_increase)
data_type <- c("H3K36m3","RNA")
names(increase_upset)<-data_type
upset(fromList(increase_upset),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      sets=data_type,nsets = 100,     # 绘制的最大集合个数
      nintersects = 20, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 1.8 # 文字标签的大小
)
decrease_upset <- list(H3K36me3_decrease,RNA_decrease)
names(decrease_upset)<-data_type
upset(fromList(decrease_upset),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      sets=data_type,nsets = 100,     # 绘制的最大集合个数
      nintersects = 20, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 1.8 # 文字标签的大小
)
