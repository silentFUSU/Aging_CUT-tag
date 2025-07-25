rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
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
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum")
resolution <- "10000"
antibody <- "H3K9me3"
condition <- "out-out"
broad_peak_min_length <- 200000
summary <- data.frame()
for(tissue in tissues){
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  if(antibody %in% c("H3K9me3","H3K36me3","H3K27me3")){
    df <- read.csv(paste0("result/HiC/",tissue,"/with_histone/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))  
  }else{
    df <- read.csv(paste0("result/HiC/",tissue,"/with_histone/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))  
  }
  proportion_summary <- df %>%
    group_by(sample, distance) %>%
    mutate(total_scale_sum = sum(scale_sum)) %>%          
    group_by(sample, distance, condition) %>%
    summarize(
      condition_scale_sum = sum(scale_sum),                 
      proportion = condition_scale_sum / total_scale_sum
    ) %>%
    ungroup()
  proportion_summary <- as.data.frame(proportion_summary)
  proportion_summary <- proportion_summary[which(proportion_summary$condition==condition),]
  proportion_summary <- merge(proportion_summary, HiC_search_table,by.x="sample",by.y="sample_name")
  
  average_age_summary <- proportion_summary %>%
    group_by(age, condition, distance) %>%
    summarise(proportion_avg = mean(proportion, na.rm = TRUE))
  
  average_age_summary <- as.data.frame(average_age_summary)
  df_24M <- average_age_summary %>% filter(age == "24M")
  df_3M <- average_age_summary %>% filter(age == "3M")
  merged_df <- df_24M %>%
    inner_join(df_3M, by = c("condition", "distance"), suffix = c("_24M", "_3M"))
  log2_ratio_df <- merged_df %>%
    mutate(log2_ratio = log2(proportion_avg_24M / proportion_avg_3M))
  log2_to_plot <- log2_ratio_df[,c("condition","distance","log2_ratio")]
  log2_to_plot <- log2_to_plot[which(log2_to_plot$distance >=100 & log2_to_plot$distance <= 6000),]
  # log2_to_plot <- log2_to_plot[which(log2_to_plot$distance >0),]
  averaged_log2_to_plot <- log2_to_plot %>%
    mutate(distance_group = floor(distance / 100) * 100) %>%  
    group_by(distance_group) %>%
    summarize(mean_log2_ratio = mean(log2_ratio, na.rm = TRUE))  
  averaged_log2_to_plot$distance_group[which(averaged_log2_to_plot$distance_group==0)] <- 1
  averaged_log2_to_plot$labels = ifelse(as.numeric(resolution) * averaged_log2_to_plot$distance_group >= 1e6, 
                                        paste0(round(as.numeric(resolution) * averaged_log2_to_plot$distance_group / 1000000, 1), "M"), 
                                        ifelse(as.numeric(resolution) * averaged_log2_to_plot$distance_group >= 1e3, 
                                               paste0(round(as.numeric(resolution) * averaged_log2_to_plot$distance_group / 1000, 1), "K"),
                                               as.numeric(resolution) * averaged_log2_to_plot$distance_group))
  averaged_log2_to_plot <- averaged_log2_to_plot[,c("labels","mean_log2_ratio")]
  colnames(averaged_log2_to_plot)[2] <- tissue_label_change(tissue)
  
  if(nrow(summary)==0){
    summary <- averaged_log2_to_plot 
  }else{
    summary <- merge(summary,averaged_log2_to_plot,by="labels")
  }
}
to_plot <- summary 
to_plot$labels <- factor(to_plot$labels,levels = c("10K",paste0(1:200,"M")))
to_plot <- to_plot[order(to_plot$labels),]
rownames(to_plot) <- to_plot$labels
to_plot <- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-0.2, -0.06, length.out = 40), seq(-0.05, 0.05, length.out = 20), seq(0.06, 0.2, length.out = 40))  
# tissue_order <- c("Lung","Cortex","Hippocampus","Cerebellum","Heart","Stomach","Liver","Thymus","Kidney","Mammary Gland","Skin","Bone Marrow","Colon")
# tissue_order <- c("Kidney","Cerebellum","Stomach","Cortex","Skin","Hippocampus","Mammary Gland","Lung","Liver","Heart","Thymus","Muscle","Bone Marrow","Colon")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,breaks = breaks, color = color_palette,clustering_distance_cols = "manhattan",main = paste0(antibody," ",condition))

if(condition !="within-within"){
  to_plot <- to_plot[,tissue_order]
  pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,breaks = breaks, color = color_palette,clustering_distance_cols = "manhattan",main = paste0(antibody," ",condition))
}else{
  p <- pheatmap::pheatmap(to_plot,cluster_rows = F,breaks = breaks, color = color_palette,clustering_distance_cols = "manhattan",main = paste0(antibody," ",condition))
  tissue_order <- p[["gtable"]][["grobs"]][[4]][["label"]]
  }

