rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(deepToolsDownstream)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)  
tissue <- args[1]  
condition <- args[2]
se <- importCount(paste0("result/WGBS/",tissue,"/WGBS_change_in_DEG/matrix/WGBS_change_in_DEG_body_",condition,".mat.gz"))
names(se@assays@data@listData) <-c("young1","young2","old1","old2")
se@metadata$sample_labels <-c("young1","young2","old1","old2")
cols  <- setNames(c("red","red","blue","blue"),c("old1","old2","young1","young2"))
p <- plotProfile(se,facet = NULL) +scale_x_continuous(  
  labels = c("-10000", "TSS", "TES", "10000")) +
  scale_color_manual(values = cols) +
  theme(text = element_text(size = 20)) +
  ggtitle(paste0(tissue,"\nCG% in gene expression ",condition," region"))

ggsave(paste0("result/WGBS/",tissue,"/WGBS_change_in_DEG/plot/WGBS_change_in_DEG_body_",condition,".png"),p,width = 8,height = 6,type="cairo")
