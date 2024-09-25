rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
# tissues <- c("lung","liver","kidney","mammarygland")
tissues <- c("ileum","Hip")
# summary <- data.frame(sample = as.character(),
#                       CG = as.numeric(),
#                       tissue = as.character())
summary <- read.csv("data/samples/WGBS/CG_manual.csv")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  t_search_table <- search_table[which(search_table$tissue == tissue),]
  samples <- t_search_table$sample_name
  for(sample in samples){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    CG <- sum(df$V4)/sum(df$V5)
    t_summary <- data.frame(sample = sample, CG = CG,tissue = tissue)
    summary <- rbind(summary,t_summary)
  }
}
write.csv(summary,"data/samples/WGBS/CG_manual.csv",row.names = F)
