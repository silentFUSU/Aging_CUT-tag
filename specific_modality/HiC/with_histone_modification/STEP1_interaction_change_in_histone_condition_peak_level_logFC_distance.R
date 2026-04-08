rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
library(ggsignif)
library(GenomicRanges)
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
tissue <- "ileum"
antibody <- "H3K9me3"
resolution <- "10000" 
broad_peak_min_length <- 200000
interaction_change_in_histone_condition_peak_level_logFC_distance <- function(tissue,resolution,antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    window_size <- "5000"
    gap_size <- "10000"
    histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"))
    histone <- histone[which((histone$V3-histone$V2 +1)>broad_peak_min_length),]
  }else{
    histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_young_old_narrowpeak.bed"))
  }
  
  histone$peaks <- paste0("peaks",c(1:nrow(histone)))
  histone$V2 <- histone$V2+1
  histone <- as.data.table(histone)
  setDT(histone)
  setkey(histone,V1,V2,V3)
  
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed$V2 <- bed$V2+1
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  overlaps <- foverlaps(histone, bed, type = "any", nomatch = 0L)  
  overlaps$condition <- "within"
  
  bed$condition <- "out"
  bed$condition[which(bed$V4 %in% overlaps$V4)] <- "within"
  
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
    df[, scale := (V3 / sum(V3)) * 1000000]
    
    df[, condition := paste(condition1, condition2, sep = "-")]
    df[condition == "within-out", condition := "out-within"]
    
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
  conditions <- c("within-within","out-within","out-out")
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
  combined_plot <- plot_a_list(p_list,no_of_cols = 1,no_of_rows = 3)
  dir.create(paste0("result/HiC/",tissue,"/with_histone/"))
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    dir.create(paste0("result/HiC/",tissue,"/with_histone/peaks_larger_",broad_peak_min_length,"/"))
    ggsave(paste0("result/HiC/",tissue,"/with_histone/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_relationship_with_",antibody,"_distance.png"),combined_plot,width=10,height = 15,type="cairo")
  }else{
    ggsave(paste0("result/HiC/",tissue,"/with_histone/",resolution,"_interaction_relationship_with_",antibody,"_distance.png"),combined_plot,width=10,height = 15,type="cairo")
  }

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
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    write.csv(summary,paste0("result/HiC/",tissue,"/with_histone/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))
    ggsave(paste0("result/HiC/",tissue,"/with_histone/peaks_larger_",broad_peak_min_length,"/",resolution,"_interaction_log2Foldchange_relationship_with_",antibody,"_distance.png"),p,width=10,height = 6,type="cairo")
  }else{
    write.csv(summary,paste0("result/HiC/",tissue,"/with_histone/",resolution,"_interaction_relationship_with_",antibody,"_distance.csv"))
    ggsave(paste0("result/HiC/",tissue,"/with_histone/",resolution,"_interaction_log2Foldchange_relationship_with_",antibody,"_distance.png"),p,width=10,height = 6,type="cairo")
  }
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle")
tissues <- c("pancreas","spleen","ileum","cecum")
antibodys <- c("H3K9me3","H3K27me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")

for(tissue in tissues){
  for(antibody in antibodys){
     interaction_change_in_histone_condition_peak_level_logFC_distance(tissue,resolution,antibody)
  }
}
