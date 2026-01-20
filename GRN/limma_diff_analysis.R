rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyr)
library(dplyr)
library(tidyverse)
library(plotly)
library(limma)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
for(tissue in tissues){
  df <- read.table("data/samples/GRN/TF_pagerank_sample_new.txt")
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue_label==tissue),]
  
  df <- df[,search_table$sample_name[which(search_table$tissue_label==tissue)]]
  # df<- df[rowSums(df) != 0, ]
  group <- factor(search_table$age, levels = c("3m", "24m"))
  design <- model.matrix(~0 + group)
  colnames(design) <- c("young", "old")
  fit <- lmFit(df, design)
  contrast <- makeContrasts(old - young, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  results <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")
  results$tissue <- tissue
  dir.create("data/samples/GRN/TF_pagerank_limma_remove_zero_row/",recursive = T)
  write.table(results,paste0("data/samples/GRN/TF_pagerank_limma_remove_zero_row/",tissue,"_TF_diff.txt"),append = F,quote = F,sep = "\t")
}
