rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(deepToolsDownstream)
library(ggplot2)
tissue <- "testis"
se <- importCount(paste0("result/",tissue,"/H3K9me3/matrix/testis_H3K9me3_in_H3K9me3-W5000-G10000-E100_filter_peaks.mat.gz"))
plotProfile(se,facet = NULL) +scale_x_continuous(  
  labels = c("-10000", "TSS", "TES", "10000")) +
  theme(text = element_text(size = 15)) +
  ylab("Delta")+
  geom_hline(yintercept = 0, linetype="dashed", color = "red")+
  ggtitle(paste0(tissue,"\nCG% Delta in gene expression ",condition," region"))
