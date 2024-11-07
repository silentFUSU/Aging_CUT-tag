rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
df <- read.table(paste0("data/samples/HiC/lung/ob_ex_matrix/WJH-100-Lung/WJH-100-Lung_50000_chr1_ob_ex_Matrix_test.txt"),sep = "\t")
colnames(df) <- df[1,]
df <- df[-1,]
rownames(df) <- df[,1]
df <- df[,-c(1:2)]
df_numeric <- as.data.frame(lapply(df, function(col) as.numeric(as.character(col))))  
pheatmap::pheatmap(df_numeric,cluster_cols = F,cluster_rows = F,show_rownames = F,show_colnames = F,na_col = "white")
df_log2 <- as.data.frame(lapply(df_numeric, function(col) {  
  transformed <- log2(col)  
  transformed[is.infinite(transformed)] <- NA  
  return(transformed)  
}))  
pheatmap::pheatmap(df_log2,cluster_cols = F,cluster_rows = F,show_rownames = F,show_colnames = F,na_col = "white",breaks = seq(-4, 4, length.out=101))
