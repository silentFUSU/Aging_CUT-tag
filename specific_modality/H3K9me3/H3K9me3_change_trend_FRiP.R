rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)
library(DiffBind)

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
tissue <- "testis"
antibody <- "H3K9me3"
window_size <- "5000"
gap_size <- "10000"
FRiP_trend <- function(tissue, antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
  samples <- search_table$sample_name
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"))
  peaks <- peaks[which((peaks$V3-peaks$V2)>100000),]
  write.table(peaks,paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_filter.bed"),quote = F,row.names = F,col.names = F,sep = "\t")
  # system(paste0("bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K9me3/H3K9me3_change_trend_FRiP.sh ",tissue),intern=T)
  
  }

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
for(tissue in tissues){
  FRiP_trend(tissue,antibody)
}
