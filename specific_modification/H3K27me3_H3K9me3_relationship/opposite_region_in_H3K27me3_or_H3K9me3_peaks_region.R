rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(ggrepel)
library(ggalluvial)  
library(scales)  
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

tissue <- "kidney"
opposite_change_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed"))
H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
H3K27me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_10kb_in_old_merge-W1000-G3000-E100.bed"))

opposite_change_region <- opposite_change_region[which(opposite_change_region$V4 %in% c(H3K9me3_peak_region$V4,H3K27me3_peak_region$V4)),]





