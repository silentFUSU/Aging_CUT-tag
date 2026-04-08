rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(corrplot)
tissues <- c("brain","CB","stomach","colon","lung","liver","thymus","heart","bonemarrow","kidney")
resolution <- "10000"
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
tissues_interaction_change_similarity <- function(tissues, resolution){
  tissues_diff <- data.frame()
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/HiC_search_table.csv")
    search_table <- search_table[which(search_table$tissue == tissue),]
    bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",search_table$sample_name[1],"_",resolution,"_abs.bed"))
    bed <- bed[which(bed$V1 %in% paste0("chr",c(1:19,"X"))),]
    young_samples <- search_table$sample_name[which(search_table$age=="3M")]
    old_samples <- search_table$sample_name[which(search_table$age=="24M")]
    young_counts <- data.frame()
    old_counts <- data.frame()
    bed$label <- paste(bed$V1,bed$V2,bed$V3,sep = "-")
    for(sample in young_samples){
      counts <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,".matrix"))
      filtered_counts <- counts[abs(V2 - V1) > 9 & abs(V2 - V1) <= 600]  
      filtered_counts[, label := paste0(V1,"-",V2)]
      filtered_counts <- filtered_counts[, .(label, V3)]  
      filtered_counts <- as.data.frame(filtered_counts)
      if(nrow(young_counts)==0){
        young_counts <- filtered_counts
      }else{
        young_counts <- merge(young_counts,filtered_counts,by="label")
      }
    }
    young_counts$row_means <- rowMeans(young_counts[,-1])
    
    for(sample in old_samples){
      counts <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,".matrix"))
      filtered_counts <- counts[abs(V2 - V1) > 9 & abs(V2 - V1) <= 600]  
      filtered_counts[, label := paste0(V1,"-",V2)]
      filtered_counts <- filtered_counts[, .(label, V3)]  
      filtered_counts <- as.data.frame(filtered_counts)
      if(nrow(old_counts)==0){
        old_counts <- filtered_counts
      }else{
        old_counts <- merge(old_counts,filtered_counts,by="label")
      }
    }
    old_counts$row_means <- rowMeans(old_counts[,-1])
    tissue_diff <- merge(young_counts,old_counts,by="label")
    tissue_diff$diff <- log2(tissue_diff$row_means.y/tissue_diff$row_means.x *(mean(tissue_diff$row_means.x)/mean(tissue_diff$row_means.y)))
    tissue_diff <- tissue_diff[,c("label","diff")]    
    colnames(tissue_diff)[2] <- tissue_label_change(tissue)
    
    if(nrow(tissues_diff)==0){
      tissues_diff <- tissue_diff
    }else{
      tissues_diff <- merge(tissues_diff,tissue_diff,by="label")
    }
  }
  return(tissues_diff)
}
tissues_diff <- tissues_interaction_change_similarity(tissues,resolution)
write.csv(tissues_diff,"data/samples/HiC/tissues_diff_similarity.csv",quote = F,row.names = F)

tissues_diff <- read.csv("data/samples/HiC/tissues_diff_similarity.csv",row.names = 1)
cor <- cor(tissues_diff)
diag(cor) <- NA  
pheatmap::pheatmap(cor)
