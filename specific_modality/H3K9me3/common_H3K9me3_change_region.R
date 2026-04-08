rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
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
common_increase <- read.csv("data/samples/all/H3K9me3/common_increase_10kb_bins_after_remove_batch_effect.csv")
common_decrease <- read.csv("data/samples/all/H3K9me3/common_decrease_10kb_bins_after_remove_batch_effect.csv")
colnames(common_increase)[2] <- "count"
colnames(common_decrease)[2] <- "count"
ggplot(common_increase, aes(x = count)) +  
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +  
  labs(x = "Count", y = "频数") +  
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.2) +
  ggtitle(paste0("H3K9me3 increase common region")) +
  theme_minimal() +
  theme(  
    axis.title.x = element_text(size = 16),  # Increase x-axis label size  
    axis.title.y = element_text(size = 16),  # Increase y-axis label size  
    axis.text = element_text(size = 14),     # Increase axis tick labels size  
    plot.title = element_text(size = 18)     # Increase plot title size  
  )

ggplot(common_decrease, aes(x = count)) +  
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +  
  labs(x = "Count", y = "频数") +  
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.2) +
  ggtitle(paste0("H3K9me3 decrease common region")) +
  theme_minimal() +
  theme(  
    axis.title.x = element_text(size = 16),  # Increase x-axis label size  
    axis.title.y = element_text(size = 16),  # Increase y-axis label size  
    axis.text = element_text(size = 14),     # Increase axis tick labels size  
    plot.title = element_text(size = 18)     # Increase plot title size  
  )
