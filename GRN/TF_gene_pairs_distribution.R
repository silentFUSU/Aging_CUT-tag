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
load("data/samples/GRN/grn_union_tissue.rdata")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
summary <- data.frame()
for(tissue in tissues){
  df <- grn_tissue[[tissue]]
  df <- as.data.frame(table(df$TF))
  df$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,df)
}
to_plot <- summary
ggplot(to_plot, aes(x = Freq, fill = tissue)) +
  geom_density(alpha = 0.2) +
  labs(x = "TF_GENE_pairs", y = "Densixty") +
  theme_minimal() + xlim(0,500)+
  facet_wrap(~ tissue, scales = "free")
mean(to_plot$Freq)
median(to_plot$Freq)

### same trend
Up_summary <- data.frame()
Down_summary <- data.frame()
for(tissue in tissues){
  df <- grn_tissue[[tissue]]
  df$tf_logFC <- log2(df$tf_old/df$tf_young)
  df$gene_logFC <- log2(df$gene_old/df$gene_young)
  df$peak_logFC <- log2(df$peak_old/df$peak_young)
  df$condition <- "other"
  df$condition[which(df$gene_logFC >0 & df$tf_logFC >0 & df$peak_logFC >0)] <- "Up"
  df$condition[which(df$gene_logFC <0 & df$tf_logFC <0 & df$peak_logFC <0)] <- "Down"
  
  Up <- as.data.frame(table(df$TF[which(df$condition=="Up")]))
  Down <- as.data.frame(table(df$TF[which(df$condition=="Down")]))
  Up$tissue <- tissue_label_change(tissue)
  Down$tissue <- tissue_label_change(tissue)
  Up_summary <- rbind(Up_summary,Up)
  Down_summary <- rbind(Down_summary,Down)
}
mean(Up_summary$Freq)
median(Up_summary$Freq)
ggplot(Up_summary, aes(x = Freq, fill = tissue)) +
  geom_density(alpha = 0.2) +
  labs(x = NULL, y = "Density", title = "Up trend") +
  theme_minimal() +
  facet_wrap(~ tissue, scales = "free")

mean(Down_summary$Freq)
median(Down_summary$Freq)
ggplot(Down_summary, aes(x = Freq, fill = tissue)) +
  geom_density(alpha = 0.2) +
  labs(x = NULL, y = "Density", title = "Down trend") +
  theme_minimal() +
  facet_wrap(~ tissue, scales = "free")

