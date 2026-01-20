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
df <- read.table("data/samples/GRN/TF_pagerank_sample_new.txt")
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
df <- df[,search_table$sample_name[which(search_table$tissue=="Lung")]]
df <- df[rowSums(df != 0) > 0, ]

df_long <- df %>%
  rownames_to_column("TF") %>%
  pivot_longer(cols = -TF, names_to = "sample_name", values_to = "value")

df_long <- merge(df_long,search_table,by="sample_name")
# df_long$value <- log2(df_long$value)
df_long$value <- log(df_long$value+0.0000001)

to_plot <- df_long
ggplot(to_plot, aes(x = value, fill = age)) +
  geom_density(alpha = 0.2) +
  labs(x = "Expression Level", y = "Densixty", title = "Density Plot of Expression Levels by Tissue") +
  theme_minimal() +
  facet_wrap(~ tissue, scales = "free")
