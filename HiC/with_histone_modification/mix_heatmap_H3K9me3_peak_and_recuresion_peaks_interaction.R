rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
library(ggsignif)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect",axes = "collect_x")
}
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

df1 <- read.csv("data/samples/HiC/all/interaction_change_in_histone_condition_peak_level_condition_proportion_logFC_distance_median_heatmap.csv")
df2 <- read.csv("data/samples/HiC/all/interaction_change_in_histone_condition_H3K9me3_between_kmeans_recursion_peak_level_logFC_distance_median_heatmap.csv")
to_plot <- merge(df1,df2,by="X")
to_plot <- to_plot[,c("X","out.out.x","out.within","within.within",
                      "kmeans1.kmeans1","kmeans2.kmeans2","Stable.Stable","kmeans3.kmeans3",
                      "kmeans4.kmeans4","out.out.y")]
colnames(to_plot) <- c("tissue","out-out","out-within","within-within",
                       "kmeans1-kmeans1","kmeans2-kmeans2","kmeans3-kmeans3",
                       "kmeans4-kmeans4","Stable-Stable","other-other")
rownames(to_plot) <- to_plot$tissue
to_plot <- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-0.2, -0.06, length.out = 40), seq(-0.05, 0.05, length.out = 20), seq(0.06, 0.2, length.out = 40)) 
to_plot <- to_plot[order(to_plot$`within-within`),]
pheatmap::pheatmap(to_plot,breaks = breaks,color = color_palette,cluster_rows = F,cluster_cols = F,border_color = "black",gaps_col = 3,filename = "result/figures/HiC_H3K9me3_peaks_interaction.pdf",height = 8,width = 6)
