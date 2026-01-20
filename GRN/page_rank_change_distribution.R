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
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "ovary"
df <- read.table("data/samples/GRN/TF_pagerank_sample_new.txt")
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
search_table <- search_table[which(search_table$tissue_label==tissue),]

df <- df[,search_table$sample_name]
df_zscore <- as.data.frame(scale(df))
limma_result <- read.table("data/samples/GRN/TF_pagerank_limma/ovary_TF_diff.txt")

### raw value distribution
young_df <- df[,search_table$sample_name[which(search_table$age=="3m")]]
young_df$young <- rowMeans(young_df)
old_df <- df[,search_table$sample_name[which(search_table$age=="24m")]]
old_df$old <- rowMeans(old_df)

to_plot <- merge(young_df[,"young",drop=F],old_df[,"old",drop=F],by="row.names")
to_plot$logFC <- log2(to_plot$old/to_plot$young)
to_plot$delta <- to_plot$old - to_plot$young

ggplot(to_plot, aes(x = delta)) +
  geom_density(alpha = 0.2) +
  labs(x = "raw data O - Y", y = "Densixty") +
  theme_minimal() 

ggplot(to_plot, aes(x = logFC)) +
  geom_density(alpha = 0.2) +
  labs(x = "raw data log2(O/Y)", y = "Densixty") +
  theme_minimal() 

### zscore distribution
young_df <- df_zscore[,search_table$sample_name[which(search_table$age=="3m")]]
young_df$young <- rowMeans(young_df)
old_df <- df_zscore[,search_table$sample_name[which(search_table$age=="24m")]]
old_df$old <- rowMeans(old_df)

to_plot <- merge(young_df[,"young",drop=F],old_df[,"old",drop=F],by="row.names")
to_plot$delta <- to_plot$old - to_plot$young

ggplot(to_plot, aes(x = delta)) +
  geom_density(alpha = 0.2) +
  labs(x = "zscore data O - Y", y = "Densixty") +
  theme_minimal() 

### limma logFC distribution
limma_result <- read.table("data/samples/GRN/TF_pagerank_limma/ovary_TF_diff.txt")

ggplot(limma_result, aes(x = logFC)) +
  geom_density(alpha = 0.2) +
  labs(x = "limma logFC", y = "Densixty") +
  theme_minimal() 



