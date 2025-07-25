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
tissue <- "cecum"
antibody <- "H3K9me3"
resolution <- "10000"
broad_peak_min_length <- 200000
interaction_change_in_histone_condition_proportion_peak_level_logFC_distance <- function(tissue,resolution,antibody){
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
  conditions <- c("within-within","out-within","out-out")
  color <- setNames(HiC_search_table$color,HiC_search_table$sample_name)
  
  breaks = seq(0,7000,1000)
  labels = ifelse(as.numeric(resolution) * abs(breaks) > 1e6, 
                  paste0(round(as.numeric(resolution) * abs(breaks) / 1000000, 1), "M"), 
                  ifelse(as.numeric(resolution) * abs(breaks) > 1e3, 
                         paste0(round(as.numeric(resolution) * abs(breaks) / 1000, 1), "K"),
                         as.numeric(resolution) * abs(breaks)))
  
  p_list <- list()
  for(condition in conditions){
    p_list[[condition]] <- ggplot(proportion_summary[which(proportion_summary$condition==condition & proportion_summary$distance >=100 & proportion_summary$distance <= 6000),], aes(x = distance, y = proportion, color =sample)) +
      geom_point() +
      geom_line() +
      scale_color_manual(values = color)+
      labs(x = "Distance",
           y = "Proportion") +
      ggtitle(paste0(tissue_label_change(tissue)," ",condition," ",antibody," peaks"))+
      scale_x_continuous("Distance", breaks=breaks, labels=labels )+
      theme_minimal()+
      theme(
        text = element_text(size = 14),  # Increase overall text size
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text
      )
  }
  combined_plot <- plot_a_list(p_list,no_of_cols = 1,no_of_rows = 3)
  dir.create(paste0("result/HiC/",tissue,"/with_histone/condition_proportion/"))
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    dir.create(paste0("result/HiC/",tissue,"/with_histone/condition_proportion/peaks_larger_",broad_peak_min_length,"/"))
    ggsave(paste0("result/HiC/",tissue,"/with_histone/condition_proportion/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_relationship_with_",antibody,"_condition_proportion_distance.png"),combined_plot,width=10,height = 15,type="cairo")
  }else{
    ggsave(paste0("result/HiC/",tissue,"/with_histone/condition_proportion/",resolution,"_interaction_relationship_with_",antibody,"_condition_proportion_distance.png"),combined_plot,width=10,height = 15,type="cairo")
  }

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
  p <- ggplot(log2_to_plot[which(log2_to_plot$distance >=100 & log2_to_plot$distance <= 6000),], aes(x = distance, y =log2_ratio, color=condition)) +
    geom_point(alpha = 0.2) +
    geom_smooth(se = FALSE, method = "loess") + 
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(
      x = "Distance",
      y = "log2(Fold change)") +
    ggtitle(paste0(tissue_label_change(tissue)," condition proportion change relationship with ",antibody," peaks"))+
    scale_x_continuous("Distance", breaks=breaks, labels=labels )+
    theme_minimal()+
    theme(
      text = element_text(size = 14),  # Increase overall text size
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text
    )
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    ggsave(paste0("result/HiC/",tissue,"/with_histone/condition_proportion/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_log2Foldchange_relationship_with_",antibody,"_condition_proportion_distance.png"),p,width=10,height = 6,type="cairo")
  }else{
    ggsave(paste0("result/HiC/",tissue,"/with_histone/condition_proportion/",resolution,"_interaction_log2Foldchange_relationship_with_",antibody,"_condition_proportion_distance.png"),p,width=10,height = 6,type="cairo")
  }
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle")
antibodys <- c("H3K27me3","H3K9me3","H3K27ac","H3K4me3","H3K4me1","H3K36me3")
for(antibody in antibodys){
  for(tissue in tissues){
    interaction_change_in_histone_condition_proportion_peak_level_logFC_distance(tissue,resolution,antibody)
  }
}
