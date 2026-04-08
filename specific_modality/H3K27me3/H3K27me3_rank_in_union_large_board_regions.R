rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
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
logFC_rank <- data.frame() 
common_increase <- read.csv("data/samples/all/H3K27me3/common_increase_10kb_bins_after_remove_batch_effect.csv")
common_decrease <- read.csv("data/samples/all/H3K27me3/common_decrease_10kb_bins_after_remove_batch_effect.csv")
regions <- c(common_increase$Geneid[which(common_increase$n>1)],common_decrease$Geneid[which(common_decrease$n>1)])

for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_10kb_in_young_old_merge-W5000-G10000-E100_bedtools_filtered.bed"))
  # df <- df[which(df$Geneid %in% peaks$V4 & df$Significant != "Stable"),]
  df <- df[which(df$Geneid %in% peaks$V4 & df$Geneid %in% regions),]
  mean_logFC <- median(df$LogFC.old.young)
  t_logFC_rank <- data.frame(tissue=tissue_label_change(tissue),mean_logFC=mean_logFC,counts=nrow(df))  
  logFC_rank <- rbind(logFC_rank,t_logFC_rank)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(logFC_rank$tissue)))
logFC_rank <- logFC_rank[order(logFC_rank$mean_logFC),]
logFC_rank$tissue <- factor(logFC_rank$tissue,levels=logFC_rank$tissue)
ggplot(logFC_rank,mapping = aes(x=mean_logFC,y=tissue,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("median log2(Fold Change)")+ggtitle(paste0("Union broad regions"))+
  theme(text = element_text(size = 18))+ scale_fill_manual(values = color) +guides(fill= guide_legend(title = ""))+
  geom_text(data = logFC_rank[which(logFC_rank$mean_logFC >0),],aes(label = counts), position = position_dodge2(width = 0.9), hjust = 1, size = 3)+
  geom_text(data = logFC_rank[which(logFC_rank$mean_logFC <0),],aes(label = counts), position = position_dodge2(width = 0.9), hjust = -0.1, size = 3)

