rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
tissues <- "ovary"
df <- read.csv("data/samples/CB/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv")
df2 <- read.csv("data/samples/CB/H3K9me3/H3K9me3_10kb_bins_diff.csv")
df <- df[which(df$Significant_bar !="Stable"),]
df2 <- df2[which(df2$Significant_bar != "Stable"),]
not_significant <- df[which(! df$Geneid %in% df2$Geneid),]
df$condition <- "not_significant"
df$condition[which(df$Geneid %in% df2$Geneid)] <- "significant"
to_plot <- df[,c("FDR.old.young","Significant_bar","condition")]
to_plot$log10fdr <- -log10(to_plot$FDR.old.young)
ggplot(to_plot, aes(x = Significant_bar, y = log10fdr, fill = condition)) +  
  geom_boxplot() +  
  theme_minimal() +  
  labs(title = "siginificant or not before removing batch effect",  
       x = NULL,  
       y = "-log10(fdr)") +  
  scale_fill_brewer(palette = "Pastel1")  

tissues <- "ovary"
df <- read.csv("data/samples/ovary/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv")
df2 <- read.csv("data/samples/ovary/H3K9me3/H3K9me3_10kb_bins_diff.csv")
df <- df[which(df$Significant_bar !="Stable"),]
df2 <- df2[which(df2$Significant_bar != "Stable"),]

df2$condition <- "not_significant"
df2$condition[which(df2$Geneid %in% df$Geneid)] <- "significant"
to_plot <- df2[,c("FDR.old.young","Significant_bar","condition")]
not_significant <- df2[which(! df2$Geneid %in% df$Geneid),]
to_plot$log10fdr <- -log10(to_plot$FDR.old.young)
ggplot(to_plot, aes(x = Significant_bar, y = log10fdr, fill = condition)) +  
  geom_boxplot() +  
  theme_minimal() +  
  labs(title = "siginificant or not before removing batch effect",  
       x = NULL,  
       y = "-log10(fdr)") +  
  scale_fill_brewer(palette = "Pastel1")  



