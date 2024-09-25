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
library(DSS)
tissues <- c("liver","lung","mammarygland","kidney")
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }
  }
  return(tissue_label)
}
search_table <-read.csv("data/samples/all/WGBS_search_table.csv")
tissue <- "lung"
bin_size <- "10kb"
WGBS_change_compress_to_bin <- function(tissue,bin_size){
  t_search_table <- search_table[which(search_table$tissue == tissue),]
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/",bin_size,"_bin_level_CpG.csv"))
  colnames(t_search_table)[3] <- "sample"
  df <- merge(df,t_search_table,by = "sample")
  switch_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant.bed"))
  df <- df[which(df$bin_id %in% switch_region$V4),]
  df$age <- factor(df$age,levels=c("3M","24M"))
  t <- t.test(df$percent[which(df$age=="3M")],df$percent[which(df$age=="24M")])
  p <- ggplot(df, aes(x = age, y = percent,fill=sample)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20)) +
    labs(title = paste0(tissue_label_change(tissue),"\nHeterochromatin Switching Region"), x = NULL, y = "CG%") +
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red")  # 在右下角标注 p 值  
  print(p)
  }
for(tissue in tissues){
  WGBS_change_compress_to_bin(tissue,bin_size)
}
