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
zscore <- function(x) {
  (x - mean(x)) / sd(x)
}
tissue <- "ovary"
TF <- "Nr5a2"
df <- read.table("data/samples/GRN/TF_pagerank_sample_new.txt")

search_table <- read.csv("data/samples/all/RNA_search_table.csv")
search_table <- search_table[which(search_table$tissue_label==tissue),]
df <- df[,search_table$sample_name]
df <- df[rowSums(df != 0) > 0, ]
rownames <- rownames(df)
df <- as.data.frame(lapply(df, zscore))
rownames(df) <- rownames
df <- df[TF,]

to_plot <- as.data.frame(t(df))
to_plot <- merge(to_plot, search_table,by.x="row.names",by.y="sample_name")
colnames(to_plot)[c(1:2)] <- c("sample_name","TF")

to_plot_ave <- to_plot %>%
  group_by(age) %>%
  summarise(mean_TF = mean(TF, na.rm = TRUE))
to_plot_ave$age <- factor(to_plot_ave$age,levels=rev(c("3m","24m")))
p <- ggplot(to_plot_ave, aes(x = mean_TF, y = age,fill=age)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = paste("Ovary ",TF),
       x = "Mean TF",
       y = "Age")
ggsave(paste0("result/Sup_figures/ovary_",TF,"_pagerank_zscore.pdf"),width = 6,height = 4)
