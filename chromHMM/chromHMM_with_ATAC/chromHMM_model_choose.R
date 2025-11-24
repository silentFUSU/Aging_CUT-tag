rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
state <- "14"
model_plot <- function(tissues,state){
  model <- read.delim(paste0("result/all/ChromHMM_with_ATAC/all_tissues/",state,"_all_tissues/emissions_",state,".txt"))
  model <- model[,-1]
  model <- model[,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC")]
  rownames(model) <- paste0("state",c(1:nrow(model)))
  color_palette <- colorRampPalette(c("white", "blue"))(50) 
  
  pheatmap::pheatmap(model,cluster_cols = F,cluster_rows = F,color = color_palette)
}
