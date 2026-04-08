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
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
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
tissue <- "lung"
antibody <- "H3K9me3"
resolution <- "10000"
histone <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_names <- strsplit(histone$label, "[:-]")
histone <- data.frame(
  V1 = sapply(split_names, "[", 1),
  V2 = sapply(split_names, "[", 2),
  V3 = sapply(split_names, "[", 3),
  cluster = histone$cluster
)
histone$peaks <- paste0("peaks",c(1:nrow(histone)))
histone$V2 <- as.numeric(histone$V2)
histone$V3 <- as.numeric(histone$V3)
histone$V2 <- histone$V2+1
histone <- as.data.table(histone)
setDT(histone)
setkey(histone,V1,V2,V3)
category_order <- c("kmeans1", "kmeans2", "kmeans3", "kmeans4","Stable","out")
category_map <- setNames(seq_along(category_order), category_order)

interaction_change_in_histone_condition_peak_level_logFC <- function(tissue,resolution){
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed$V2 <- bed$V2+1
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  overlaps <- foverlaps(histone, bed, type = "any", nomatch = 0L)  
  overlaps <- as.data.frame(overlaps)
  overlaps$overlap_length <- pmin(overlaps$V3, overlaps$i.V3) - pmax(overlaps$V2, overlaps$i.V2)
  overlaps <- overlaps %>%
    group_by(V4) %>%           
    slice_max(overlap_length)  
  bed$condition <- "out"
  bed$condition[which(bed$V4 %in% overlaps$V4[which(overlaps$cluster=="kmeans1")])] <- "kmeans1"
  bed$condition[which(bed$V4 %in% overlaps$V4[which(overlaps$cluster=="kmeans2")])] <- "kmeans2"
  bed$condition[which(bed$V4 %in% overlaps$V4[which(overlaps$cluster=="kmeans3")])] <- "kmeans3"
  bed$condition[which(bed$V4 %in% overlaps$V4[which(overlaps$cluster=="kmeans4")])] <- "kmeans4"
  bed$condition[which(bed$V4 %in% overlaps$V4[which(overlaps$cluster=="Stable")])] <- "Stable"

  
  summary <- data.frame()
  for(sample in HiC_search_table$sample_name){
    df <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,".matrix"))
    setkey(bed, V4)
    setkey(df, V1)
    df[, Chr1 := bed[V1, V1]]
    df[, condition1 := bed[V1, condition]]
    
    setkey(df, V2)
    df[, Chr2 := bed[V2, V1]]
    df[, condition2 := bed[V2, condition]]
    df <- df[Chr1 == Chr2]
    unique_conditions <- unique(c(df$condition1, df$condition2))
    precomputed_conditions <- as.data.table(expand.grid(unique_conditions, unique_conditions))
    precomputed_conditions[, combined := {
      idx <- category_map[c(Var1, Var2)]
      sorted_conditions <- c(Var1, Var2)[order(idx,decreasing = T)]
      paste(sorted_conditions, collapse = "-")
    }, by = .(Var1, Var2)]

    setkey(precomputed_conditions, Var1, Var2)
    df[, condition := precomputed_conditions[J(condition1, condition2), combined]]
  
    df[, scale := (V3 / sum(V3)) * 1000000]
    df <- df[, .(V1,V2,scale,condition)]
    df[, distance := V2 - V1]
    
    result <- df[, .(scale_sum = sum(scale)), by = .(condition, distance)]
    result <- as.data.frame(result)
    result$sample <- sample
    summary <- rbind(summary,result)
  }
  summary <- merge(summary, HiC_search_table,by.x="sample",by.y="sample_name")
  
  HiC_search_table$age <- factor(HiC_search_table$age,levels=c("3M","24M"))
  HiC_search_table <- HiC_search_table[order(HiC_search_table$age),]
  HiC_search_table$color <- c("#f38181", "#ff2e63","#112d4e", "#3f72af")
  color <- setNames(HiC_search_table$color,HiC_search_table$sample_name)
  breaks = seq(0,7000,1000)
  labels = ifelse(as.numeric(resolution) * abs(breaks) > 1e6, 
                  paste0(round(as.numeric(resolution) * abs(breaks) / 1000000, 1), "M"), 
                  ifelse(as.numeric(resolution) * abs(breaks) > 1e3, 
                         paste0(round(as.numeric(resolution) * abs(breaks) / 1000, 1), "K"),
                         as.numeric(resolution) * abs(breaks)))
  conditions <- sort(unique(summary$condition))
  p_list <- list()
  for(condition in conditions){
    p_list[[condition]] <- ggplot(summary[which(summary$condition==condition & summary$distance >=100 & summary$distance <= 6000),], aes(x = distance, y = log2(scale_sum), color =sample)) +
      geom_point() +
      geom_line() +
      scale_color_manual(values = color)+
      labs(title = "Scatter plot of Scale vs Distance",
           x = "Distance",
           y = "log2(Interaction)") +
      ggtitle(paste0(tissue_label_change(tissue)," ",condition," ",antibody," peaks"))+
      scale_x_continuous("Distance", breaks=breaks, labels=labels )+
      theme_minimal()+
      theme(
        text = element_text(size = 14),  # Increase overall text size
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text
      )
  }
  combined_plot <- plot_a_list(p_list,no_of_cols = 3,no_of_rows = 7)
  dir.create(paste0("result/HiC/",tissue,"/with_histone/"))
  dir.create(paste0("result/HiC/",tissue,"/with_histone/H3K9me3_recursion_peaks/"))
  ggsave(paste0("result/HiC/",tissue,"/with_histone/H3K9me3_recursion_peaks/",resolution,"_interaction_relationship_with_",antibody,"_recursion_peaks_distance.png"),combined_plot,width=30,height = 35,type="cairo")
  
  average_age_summary <- summary %>%
    group_by(age, condition, distance) %>%
    summarise(scale_avg = mean(scale_sum, na.rm = TRUE))
  average_age_summary <- as.data.frame(average_age_summary)
  df_24M <- average_age_summary %>% filter(age == "24M")
  df_3M <- average_age_summary %>% filter(age == "3M")
  merged_df <- df_24M %>%
    inner_join(df_3M, by = c("condition", "distance"), suffix = c("_24M", "_3M"))
  log2_ratio_df <- merged_df %>%
    mutate(log2_ratio = log2(scale_avg_24M / scale_avg_3M))
  log2_to_plot <- log2_ratio_df[,c("condition","distance","log2_ratio")]
  p <- ggplot(log2_to_plot[which(log2_to_plot$distance >=100 & log2_to_plot$distance <= 6000),], aes(x = distance, y =log2_ratio, color=condition)) +
    geom_point(alpha = 0.2) +
    geom_smooth(se = FALSE, method = "loess") + 
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(
      x = "Distance",
      y = "log2(Fold change)") +
    ggtitle(paste0(tissue_label_change(tissue)," interaction change relationship with ",antibody," peaks"))+
    scale_x_continuous("Distance", breaks=breaks, labels=labels )+
    theme_minimal()+
    theme(
      text = element_text(size = 14),  # Increase overall text size
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text
    )
  write.csv(summary,paste0("result/HiC/",tissue,"/with_histone/H3K9me3_recursion_peaks/",resolution,"_interaction_relationship_with_",antibody,"_recursion_peaks_distance.csv"))
}

tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","ileum","cecum")
tissues <- c("pancreas","spleen")
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- interaction_change_in_histone_condition_peak_level_logFC(tissue,resolution)
}
combined_plot <- plot_a_list(p_list,no_of_rows = 3, no_of_cols = 4)







