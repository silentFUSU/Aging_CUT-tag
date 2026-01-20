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
df <- grn_tissue[["ovary"]]
df$tf_logFC <- log2(df$tf_old/df$tf_young)
df$gene_logFC <- log2(df$gene_old/df$gene_young)
df$peak_logFC <- log2(df$peak_old/df$peak_young)
df$condition <- "other"
df$condition[which(df$gene_logFC >0 & df$tf_logFC >0 & df$peak_logFC >0)] <- "Up"
df$condition[which(df$gene_logFC <0 & df$tf_logFC <0 & df$peak_logFC <0)] <- "Down"

Nr5a2 <- df[which(df$TF=="Nr5a2"),]
nrow(Nr5a2[which(Nr5a2$tf_logFC < 0 & Nr5a2$gene_logFC <0),])

gene <- read.csv("data/samples/RNA/ovary/diff_expression_gene_change_filter_bar.csv")
gene <- gene[which(gene$Significant=="Down"),]
nrow(Nr5a2[which(Nr5a2$tf_logFC < 0 & Nr5a2$gene_logFC <0 & Nr5a2$gene %in% gene$X),])
