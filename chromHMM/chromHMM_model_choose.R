rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
model <- read.delim("result/all/ChromHMM/until_ovary/16_until_ovary/emissions_16.txt")
model <- model[,-1]
model <- model[,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")]
rownames(model) <- paste0("state",c(1:nrow(model)))
color_palette <- colorRampPalette(c("white", "blue"))(50) 
pheatmap::pheatmap(model,cluster_cols = F,cluster_rows = F,color = color_palette)














