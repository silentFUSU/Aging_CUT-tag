rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
kmeans <- paste0("kmeans",c(1:4))
summary <- data.frame()
motif <- c()
for(kmean in kmeans){
  df <- read.delim(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/motif_bg/",kmean,"/knownResults.txt"))
  df <- df[which(df$P.value < 0.05),]
  motif <- c(motif,df$Motif.Name)
  }
for(kmean in kmeans){
  df <- read.delim(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/motif_bg/",kmean,"/knownResults.txt"))
  df <- df[which(df$Motif.Name %in% motif),c("Motif.Name","P.value")]
  df$P.value <- -log10(df$P.value)
  colnames(df)[2] <- kmean
  if(nrow(summary)==0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="Motif.Name")
  }
} 
to_plot <- unique(summary)
rownames(to_plot) <- to_plot$Motif.Name
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(0, 1.2, length.out = 40), seq(1.2+0.01,1.4,length.out = 20), seq(1.4+0.01,2.6,length.out = 40))  
pheatmap::pheatmap(to_plot[,-1],breaks = breaks,color = color_palette,cluster_cols = F)

summary <- data.frame()
motif <- c()
for(kmean in kmeans){
  df <- read.delim(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/motif_without_bg/",kmean,"/knownResults.txt"))
  df <- df[which(df$P.value < 0.05),]
  motif <- c(motif,df$Motif.Name)
}
for(kmean in kmeans){
  df <- read.delim(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/motif_without_bg/",kmean,"/knownResults.txt"))
  df <- df[which(df$Motif.Name %in% motif),c("Motif.Name","P.value")]
  df$P.value <- -log10(df$P.value)
  colnames(df)[2] <- kmean
  if(nrow(summary)==0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="Motif.Name")
  }
} 
to_plot <- unique(summary)
rownames(to_plot) <- to_plot$Motif.Name
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(0, 1.2, length.out = 40), seq(1.2+0.01,1.4,length.out = 20), seq(1.4+0.01,2.6,length.out = 40))  
pheatmap::pheatmap(to_plot[,-1],breaks = breaks,color = color_palette,cluster_cols = F)
