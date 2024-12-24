rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)
options(scipen = 999)  
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
      tissue_label <- "Mammarygland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
state_num <- "14"
#scale by row
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:state_num))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:state_num))))
to_plot <- data.frame()
for(tissue in tissues){
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue_label_change(tissue)),]
  t_transfer_matrix$condition <- paste0(t_transfer_matrix$young_state,"-",t_transfer_matrix$old_state)
  t_transfer_matrix <- t_transfer_matrix %>%
    group_by(young_state) %>%
    mutate(percent_freq = Freq/sum(Freq)*100)
  t_to_plot <- data.frame(condition = t_transfer_matrix$condition, percent=t_transfer_matrix$percent_freq)
  colnames(t_to_plot)[2] <- tissue_label_change(tissue)
  if(nrow(to_plot) == 0){
    to_plot <- t_to_plot
  }else{
    to_plot <- merge(to_plot, t_to_plot, by="condition")
  }
}
rownames(to_plot) <- to_plot$condition
to_plot <- to_plot[,-1]
to_plot_cor <- cor(to_plot)
pheatmap::pheatmap(to_plot_cor,main = "Scale by row")

#scale by row remove unchange
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:state_num))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:state_num))))
to_plot <- data.frame()
for(tissue in tissues){
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue_label_change(tissue)),]
  t_transfer_matrix <- t_transfer_matrix[which(t_transfer_matrix$young_state!=t_transfer_matrix$old_state),]
  t_transfer_matrix$condition <- paste0(t_transfer_matrix$young_state,"-",t_transfer_matrix$old_state)
  t_transfer_matrix <- t_transfer_matrix %>%
    group_by(young_state) %>%
    mutate(percent_freq = Freq/sum(Freq)*100)
  t_to_plot <- data.frame(condition = t_transfer_matrix$condition, percent=t_transfer_matrix$percent_freq)
  colnames(t_to_plot)[2] <- tissue_label_change(tissue)
  if(nrow(to_plot) == 0){
    to_plot <- t_to_plot
  }else{
    to_plot <- merge(to_plot, t_to_plot, by="condition")
  }
}
rownames(to_plot) <- to_plot$condition
to_plot <- to_plot[,-1]
to_plot_cor <- cor(to_plot)
pheatmap::pheatmap(to_plot_cor,main = "Scale by row remove unchanged")
pheatmap::pheatmap(to_plot,scale = "row",main = "Scale by all condition remove empty state")

#scale by all transfer contidion
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:state_num))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:state_num))))
to_plot <- data.frame()
for(tissue in tissues){
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue_label_change(tissue)),]
  t_transfer_matrix$condition <- paste0(t_transfer_matrix$young_state,"-",t_transfer_matrix$old_state)
  t_transfer_matrix$percent_freq <- t_transfer_matrix$Freq/sum(t_transfer_matrix$Freq)*100
  t_to_plot <- data.frame(condition = t_transfer_matrix$condition, percent=t_transfer_matrix$percent_freq)
  colnames(t_to_plot)[2] <- tissue_label_change(tissue)
  if(nrow(to_plot) == 0){
    to_plot <- t_to_plot
  }else{
    to_plot <- merge(to_plot, t_to_plot, by="condition")
  }
}
rownames(to_plot) <- to_plot$condition
to_plot <- to_plot[,-1]
to_plot_cor <- cor(to_plot,method = "spearman")
pheatmap::pheatmap(to_plot_cor,main = "Scale by all condition")
pheatmap::pheatmap(to_plot,scale = "row",main = "Scale by all condition")

#scale by all transfer condition remove empty state 
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:state_num))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:state_num))))
to_plot <- data.frame()
for(tissue in tissues){
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue_label_change(tissue)),]
  t_transfer_matrix <- t_transfer_matrix[which(t_transfer_matrix$young_state!="E3" & t_transfer_matrix$young_state!="E12"),]
  t_transfer_matrix <- t_transfer_matrix[which(t_transfer_matrix$old_state!="E3" & t_transfer_matrix$old_state!="E12"),]
  t_transfer_matrix$condition <- paste0(t_transfer_matrix$young_state,"-",t_transfer_matrix$old_state)
  t_transfer_matrix$percent_freq <- t_transfer_matrix$Freq/sum(t_transfer_matrix$Freq)*100
  t_to_plot <- data.frame(condition = t_transfer_matrix$condition, percent=t_transfer_matrix$percent_freq)
  colnames(t_to_plot)[2] <- tissue_label_change(tissue)
  if(nrow(to_plot) == 0){
    to_plot <- t_to_plot
  }else{
    to_plot <- merge(to_plot, t_to_plot, by="condition")
  }
}
rownames(to_plot) <- to_plot$condition
to_plot <- to_plot[,-1]
to_plot_cor <- cor(to_plot,method = "spearman")
pheatmap::pheatmap(to_plot_cor,main = "Scale by all condition remove empty state")
color_tissues <- sapply(sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                               "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")), tissue_label_change) 
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,color_tissues)
color <- color[!is.na(names(color))]
color <- list(tissue=color)
annotation <- data.frame(tissue = colnames(to_plot))
rownames(annotation) <- annotation$tissue
pheatmap::pheatmap(to_plot,scale = "row",annotation_col = annotation,annotation_colors = color,main = "Scale by all condition remove empty state")
  