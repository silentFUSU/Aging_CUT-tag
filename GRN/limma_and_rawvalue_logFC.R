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
  pagerank_limma <- read.table(paste0("data/samples/GRN/TF_pagerank_limma_remove_zero_row_log/",tissue,"_TF_diff.txt"))
  
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue_label==tissue),]
  pagerank <- read.table("data/samples/GRN/TF_pagerank_sample_new.txt")
  pagerank <- pagerank[,search_table$sample_name]
  young_pagerank <- pagerank[,search_table$sample_name[which(search_table$age=="3m")]]
  young_pagerank$young <- rowMeans(young_pagerank)
  
  old_pagerank <- pagerank[,search_table$sample_name[which(search_table$age=="24m")]]
  old_pagerank$old <- rowMeans(old_pagerank)
  pagerank <- merge(young_pagerank[,"young",drop=F],old_pagerank[,"old",drop=F],by="row.names")
  rownames(pagerank) <- pagerank$Row.names
  pagerank <- pagerank[,-1]
  pagerank <- pagerank[rowSums(pagerank) != 0, ]
  pagerank$logFC <- log2(pagerank$old/pagerank$young)
  
  t_summary <- merge(pagerank_limma[,"logFC",drop=F],pagerank[,"logFC",drop=F],by="row.names")
  colnames(t_summary) <- c("TF","limma","raw_value")
  t_summary$TF <- paste0(tissue_label_change(tissue),"-",t_summary$TF)
  summary <- rbind(summary,t_summary)
}
to_plot <- summary
ggplot(to_plot, aes(x = limma, y = raw_value)) +
  geom_point(color = "black",alpha=0.5) +
  labs(x = "limma logFC", y = "raw_value logFC") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  )+
  xlim(-5,5)+
  ylim(-5,5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")
