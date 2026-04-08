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
  HiC_search_table$rep <- c("rep1","rep2","rep1","rep2")
  HiC_search_table$color <- c("#f38181", "#ff2e63","#112d4e", "#3f72af")
  if(antibody %in% c("H3K9me3","H3K36me3","H3K27me3")){
    summary <- read.csv(paste0("result/HiC/",tissue,"/with_histone/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))  
  }else{
    summary <- read.csv(paste0("result/HiC/",tissue,"/with_histone/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))  
  }
  reps <- c("rep1","rep2")
  log2_to_plot_summary <- data.frame()
  for(rep in reps){
    samples <- HiC_search_table$sample_name[which(HiC_search_table$rep==rep)]
    t_summary <- summary[which(summary$sample %in% samples),]
    proportion_summary <- t_summary %>%
      group_by(sample, distance) %>%
      mutate(total_scale_sum = sum(scale_sum)) %>%          
      group_by(sample, distance, condition) %>%
      summarize(
        condition_scale_sum = sum(scale_sum),                 
        proportion = condition_scale_sum / total_scale_sum
      ) %>%
      ungroup()
    proportion_summary <- merge(proportion_summary, HiC_search_table[,- which(colnames(HiC_search_table)=="color")],by.x="sample",by.y="sample_name")
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
    log2_to_plot$rep <- rep
    log2_to_plot_summary <- rbind(log2_to_plot_summary,log2_to_plot)
  }
  log2_to_plot_summary$condition <- paste0(log2_to_plot_summary$condition,"-",log2_to_plot_summary$rep)
  color <- setNames(c("#ff2e63","#fc5185","#1fab89","#62d2a2","#6eb6ff","#7098da"),c("out-out-rep1","out-out-rep2","out-within-rep1","out-within-rep2","within-within-rep1","within-within-rep2"))
  
  breaks = seq(0,7000,1000)
  labels = ifelse(as.numeric(resolution) * abs(breaks) > 1e6, 
                  paste0(round(as.numeric(resolution) * abs(breaks) / 1000000, 1), "M"), 
                  ifelse(as.numeric(resolution) * abs(breaks) > 1e3, 
                         paste0(round(as.numeric(resolution) * abs(breaks) / 1000, 1), "K"),
                         as.numeric(resolution) * abs(breaks)))
  p <- ggplot(log2_to_plot_summary[which(log2_to_plot_summary$distance >=100 & log2_to_plot_summary$distance <= 6000),], aes(x = distance, y =log2_ratio, color=condition)) +
    geom_point(alpha = 0.2) +
    geom_smooth(se = FALSE, method = "loess") + 
    scale_color_manual(values = color)+
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
  dir.create(paste0("result/HiC/",tissue,"/with_histone/condition_proportion_pseudo_replicates/"))
  ggsave(paste0("result/HiC/",tissue,"/with_histone/condition_proportion_pseudo_replicates/",resolution,"_interaction_log2Foldchange_relationship_with_",antibody,"_condition_proportion_pseudo_replicates_distance.png"),p,width=10,height = 6,type="cairo")
  }
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","cecum","ileum","pancreas","spleen","muscle","skin")
antibodys <- c("H3K27me3","H3K9me3","H3K27ac","H3K4me3","H3K4me1")
for(antibody in antibodys){
  for(tissue in tissues){
    interaction_change_in_histone_condition_proportion_peak_level_logFC_distance(tissue,resolution,antibody)
  }
}











