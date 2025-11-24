rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(patchwork)
library(edgeR)
library(corrplot) 
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
antibody <- "H3K27me3"
method <-"spearman"
options(bitmapType = "cairo")  
search_table_cut <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
search_table_atac <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
search_table <- rbind(search_table_cut,search_table_atac)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
correlation_clustering <- function(antibody,method){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  
  young_tissue_summary <- data.frame()
  old_tissue_summary <- data.frame()
  for(tissue in tissues){
    df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),header = T)
    rownames(df) <- df$Geneid
    counts <- df[,c(7:ncol(df))]
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
    t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
    CPM <- as.data.frame(edgeR::cpm(counts))
    young_CPM <- CPM[,t_search_table$sample_name[which(t_search_table$age=="3m")]]
    young_CPM$mean <- rowMeans(young_CPM)
    colnames(young_CPM)[which(colnames(young_CPM)=="mean")] <- tissue_label_change(tissue)
    young_CPM$label <- rownames(young_CPM)
    young_CPM <- young_CPM[,c("label",tissue_label_change(tissue))]
    
    old_CPM <- CPM[,t_search_table$sample_name[which(t_search_table$age=="24m")]]
    old_CPM$mean <- rowMeans(old_CPM)
    colnames(old_CPM)[which(colnames(old_CPM)=="mean")] <- tissue_label_change(tissue)
    old_CPM$label <- rownames(old_CPM)
    old_CPM <- old_CPM[,c("label",tissue_label_change(tissue))]
    
    if(nrow(young_tissue_summary)==0){
      young_tissue_summary <- young_CPM
    }else{
      young_tissue_summary <- merge(young_tissue_summary,young_CPM,by="label")
    }
    if(nrow(old_tissue_summary)==0){
      old_tissue_summary <- old_CPM
    }else{
      old_tissue_summary <- merge(old_tissue_summary,old_CPM,by="label")
    }
  }
  young_to_plot <- young_tissue_summary
  rownames(young_to_plot) <- young_to_plot$label
  young_to_plot <- young_to_plot[,-1]
  young_to_plot <- cor(young_to_plot,method = method)
  breaks <- c(seq(0.5, 0.7, length.out = 40), seq(0.71, 0.85, length.out = 20), seq(0.86, 1, length.out = 40))
  pheatmap::pheatmap(young_to_plot, breaks = breaks)
  
  old_to_plot <- old_tissue_summary
  rownames(old_to_plot) <- old_to_plot$label
  old_to_plot <- old_to_plot[,-1]
  old_to_plot <- cor(old_to_plot,method=method)
  breaks <- c(seq(0.5, 0.7, length.out = 40), seq(0.71, 0.85, length.out = 20), seq(0.86, 1, length.out = 40))
  pheatmap::pheatmap(old_to_plot, breaks = breaks)
  }

