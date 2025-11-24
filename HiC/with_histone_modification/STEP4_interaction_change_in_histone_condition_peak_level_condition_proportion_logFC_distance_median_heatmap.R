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
tissue <- "kidney"
antibody <- "H3K4me3"
resolution <- "10000"
if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
  broad_peak_min_length <- 200000
}else{
  broad_peak_min_length <- 0
}

tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum","ileum","pancreas","spleen")
median_log2FC_summary <- data.frame()
for(tissue in tissues){
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  HiC_search_table$age <- factor(HiC_search_table$age,levels=c("3M","24M"))
  HiC_search_table <- HiC_search_table[order(HiC_search_table$age),]
  HiC_search_table$color <- c("#f38181", "#ff2e63","#112d4e", "#3f72af")
  
  if(antibody %in% c("H3K9me3","H3K36me3","H3K27me3")){
    summary <- read.csv(paste0("result/HiC/",tissue,"/with_histone/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))  
  }else{
    summary <- read.csv(paste0("result/HiC/",tissue,"/with_histone/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))  
  }
  proportion_summary <- summary %>%
    group_by(sample, distance) %>%
    mutate(total_scale_sum = sum(scale_sum)) %>%          
    group_by(sample, distance, condition) %>%
    summarize(
      condition_scale_sum = sum(scale_sum),                 
      proportion = condition_scale_sum / total_scale_sum
    ) %>%
    ungroup()
  proportion_summary <- as.data.frame(proportion_summary)
  proportion_summary <- merge(proportion_summary, HiC_search_table[,-ncol(HiC_search_table)],by.x="sample",by.y="sample_name")
  average_age_proportion_summary <- proportion_summary %>%
    group_by(age, condition, distance) %>%
    summarise(proportion_avg = mean(proportion, na.rm = TRUE))
  average_age_proportion_summary <- as.data.frame(average_age_proportion_summary)
  df_24M <- average_age_proportion_summary %>% filter(age == "24M")
  df_3M <- average_age_proportion_summary %>% filter(age == "3M")
  merged_df <- df_24M %>%
    inner_join(df_3M, by = c("condition", "distance"), suffix = c("_24M", "_3M"))
  log2_ratio_df <- merged_df %>%
    mutate(log2_ratio = log2(proportion_avg_24M / proportion_avg_3M))
  log2_to_plot <- log2_ratio_df[,c("condition","distance","log2_ratio")]
  log2_to_plot <- log2_ratio_df[which(log2_ratio_df$distance >= 10000000/as.numeric(resolution) & log2_ratio_df$distance <= 40000000/as.numeric(resolution)),]
  median_log2_ratio_by_condition <- log2_to_plot %>%
    group_by(condition) %>%
    summarise(median_log2_ratio = median(log2_ratio, na.rm = TRUE))
  colnames(median_log2_ratio_by_condition)[2] <- tissue_label_change(tissue)
  
  if(nrow(median_log2FC_summary)==0){
    median_log2FC_summary <- median_log2_ratio_by_condition
  }else{
    median_log2FC_summary <- merge(median_log2FC_summary,median_log2_ratio_by_condition,by="condition")
  }
}
rownames(median_log2FC_summary) <- median_log2FC_summary$condition
median_log2FC_summary <- median_log2FC_summary[,-1]
to_plot <- as.data.frame(t(median_log2FC_summary))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-0.2, -0.06, length.out = 40), seq(-0.05, 0.05, length.out = 20), seq(0.06, 0.2, length.out = 40)) 
if(antibody=="H3K9me3"){
  p <- pheatmap::pheatmap(to_plot,breaks = breaks, border_color = "black",color = color_palette,cluster_cols = F,main = antibody,filename = "result/figures/HiC_H3K9me3_peaks_interaction.pdf",width = 6,height = 8 )
  tissue_order <- p[["gtable"]][["grobs"]][[5]][["label"]]
}else{
  to_plot <- to_plot[tissue_order,]
  pheatmap::pheatmap(to_plot,breaks = breaks,border_color = "black", color = color_palette,cluster_cols = F,filename = paste0("result/Sup_figures/HiC_",antibody,"_peaks_interaction.pdf"),width = 6,height = 8,cluster_rows = F,main = antibody)
}

