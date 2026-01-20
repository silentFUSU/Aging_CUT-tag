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
library(plotly)
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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
summary <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/GRN/TF_pagerank_limma_remove_zero_row_log/",tissue,"_TF_diff.txt"))
  df <- df[which(df$P.Value < 0.05),]  
  df$condition <- "Up"
  df$condition[which(df$logFC < 0)] <- "Down"
  t_summary <- data.frame(table(df$condition))
  t_summary$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,t_summary)
}
summary$Var1 <- factor(summary$Var1,levels = c("Up","Down"))
ggplot(summary, aes(x = tissue, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Tissue", y = "Frequency", fill = "Var1") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # 如果有太多 tissue 列，可以旋转标签
  )+
  ylim(0,200)

summary <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/GRN/TF_pagerank_limma_remove_zero_row_log/",tissue,"_TF_diff.txt"))
  df <- df[which(df$adj.P.Val < 0.05),]  
  if(nrow(df) > 0){
    df$condition <- "Up"
    df$condition[which(df$logFC < 0)] <- "Down"
    t_summary <- data.frame(table(df$condition))
    t_summary$tissue <- tissue_label_change(tissue)
  }else{
    t_summary <- data.frame(Var1=c("Down","Up"),Freq=c(0,0),tissue=tissue_label_change(tissue))
  }
  summary <- rbind(summary,t_summary)
}
summary$Var1 <- factor(summary$Var1,levels = c("Up","Down"))
p <- ggplot(summary, aes(x = tissue, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Tissue", y = "Frequency", fill = "Var1") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)  # 如果有太多 tissue 列，可以旋转标签
  )+
  ylim(0,200)
ggsave("result/Sup_figures/fdr_sig_limma_pagerank_TF.pdf",p,width = 8,height = 6)
