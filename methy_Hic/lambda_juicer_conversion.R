rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
args <- commandArgs(trailingOnly = TRUE)  
file_path <- args[1]
file_path <- "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/aligned/merged_dedup_sort_CpG.bedGraph"
data <- read.delim(file_path,skip = 1,header = F)
chrL <- data[-which(data$V1=="chrL"),]
conversion <- 100-mean(chrL$V4)
