rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
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
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
antibodys <- c("H3K27ac","H3K4me3","H3K4me1","H3K27me3","H3K36me3","ATAC","RNA")
antibody <- "H3K27ac"

for(antibody in antibodys){
  summary <- data.frame()
  for(tissue in tissues){
    if(antibody == "ATAC"){
      active_mark <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_diff_in_H3K9me3_peaks_remove_batch_effect.csv"),row.names = 1)
      active_mark <- active_mark[,c("Geneid","LogFC.old.young","logCPM","Significant")]
    }else if(antibody=="RNA"){
      active_mark <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_in_H3K9me3_peaks.csv"))
      active_mark <- active_mark[,c("X","logFC","logCPM","Significant")]
      colnames(active_mark) <-  c("Geneid","LogFC.old.young","logCPM","Significant")
    }else{
      active_mark <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_H3K9me3_peaks_remove_batch_effect.csv"),row.names = 1)
      active_mark <- active_mark[,c("Geneid","LogFC.old.young","logCPM","Significant")]
    }
    colnames(active_mark)[4] <- "active_mark_significant"
    histone <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
    histone <- histone[,c("Geneid","Length","LogFC.old.young","Significant")]
    colnames(histone)[4] <- "histone_significant"
    to_plot <- merge(histone,active_mark,by="Geneid")
    colnames(to_plot)[c(3,5)] <- c("LogFC.old.young","logFC")
    to_plot$tissue <- tissue_label_change(tissue)
    summary <- rbind(summary,to_plot)
  }
  summary <- summary[which(summary$Length > 200000),]
  to_plot <- summary %>%  
    group_by(tissue, histone_significant) %>%  
    summarise(  
      median_logFC = median(logFC, na.rm = TRUE),  
      count = n()  
    ) %>%   
    mutate(median_logFC = ifelse(count < 10, NA, median_logFC))
  to_plot <- to_plot[,-ncol(to_plot)]
  to_plot <- as.data.frame(to_plot)
  to_plot <- reshape2::dcast(to_plot, tissue ~ histone_significant, value.var = "median_logFC")  
  rownames(to_plot) <- to_plot$tissue
  to_plot <- to_plot[,-1]
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
  
  to_plot <- to_plot[,c("Up","Stable","Down")]
  
  p_value_summary <- data.frame(
    Up = rep(NA, 27), 
    Stable = rep(NA, 27), 
    Down = rep(NA, 27)  
  )
  
  rownames(p_value_summary) <- sapply(tissues, tissue_label_change)
  for(tissue in tissues){
    t_Up_summary <-  summary[which(summary$tissue == tissue_label_change(tissue) & summary$histone_significant=="Up"),]
    t_Down_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$histone_significant=="Down"),]
    t_Stable_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$histone_significant=="Stable"),]
    
    if(antibody %in% c("H3K27me3","H3K27ac","H3K4me1","ATAC","RNA")){
      if(nrow(t_Up_summary) > 10){
        test <- wilcox.test(t_Up_summary$logFC,t_Stable_summary$logFC,alternative = "less")
        p_value_summary[tissue_label_change(tissue),"Up"] <- test$p.value
      }
      if(nrow(t_Down_summary) > 10){
        test <- wilcox.test(t_Down_summary$logFC,t_Stable_summary$logFC,alternative = "greater")
        p_value_summary[tissue_label_change(tissue),"Down"] <- test$p.value
      }
    }else{
      if(nrow(t_Up_summary) > 10){
        test <- wilcox.test(t_Up_summary$logFC,t_Stable_summary$logFC,alternative = "greater")
        p_value_summary[tissue_label_change(tissue),"Up"] <- test$p.value
      }
      if(nrow(t_Down_summary) > 10){
        test <- wilcox.test(t_Down_summary$logFC,t_Stable_summary$logFC,alternative = "less")
        p_value_summary[tissue_label_change(tissue),"Down"] <- test$p.value
      }
    }
  }
  mark_significance <- function(p_value) {
    if (is.na(p_value)) {
      return(NA)
    } else if (p_value < 0.001) {
      return("***")
    } else if (p_value < 0.01) {
      return("**")
    } else if (p_value < 0.05) {
      return("*")
    } else {
      return(NA)
    }
  }
  p_value_summary <- p_value_summary %>%
    mutate(
      Up = sapply(Up, mark_significance),
      Stable = sapply(Stable, mark_significance),
      Down = sapply(Down, mark_significance)
    )
  to_plot$tissue <- rownames(to_plot)
  df_long <- reshape2::melt(to_plot)
  names(df_long) <- c("Tissue", "Type", "Value")
  
  p_value_summary$tissue <- rownames(p_value_summary)
  p_value_long <- reshape2::melt(p_value_summary,id.vars = "tissue")
  names(p_value_long) <- c("Tissue", "Type", "Label")
  merged_data <- merge(df_long, p_value_long, by = c("Tissue", "Type"), all.x = TRUE)
  
  merged_data$Value[which(merged_data$Value > 1)] <- 1
  merged_data$Value[which(merged_data$Value < -1)] <- -1
  tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary",
                    "Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
  merged_data$Tissue <- factor(merged_data$Tissue,levels = tissue_order)
  p <- ggplot(merged_data[which(merged_data$Type %in% c("Up","Down")),], aes(x = Type, y = Tissue, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-1, 1), midpoint = 0) +
    theme_minimal() +
    xlab("H3K9me3 change condition")+
    ggtitle(paste0(antibody))+
    geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0("result/Sup_figures/H3K9me3_peak_",antibody,"_change_heatmap.pdf"),p,width = 6,height = 8)
}
