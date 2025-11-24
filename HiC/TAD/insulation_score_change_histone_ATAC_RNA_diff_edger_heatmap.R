rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
tissue <- "lung"
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus","skin","muscle","cecum","ileum","pancreas","spleen")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC","RNA")
for(antibody in antibodys){
  summary <- data.frame()
  for(tissue in tissues){
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")){
      active_mark <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_20000_redundant_tads.csv"),row.names = 1)
      active_mark <- active_mark[,c("Geneid","LogFC.old.young","logCPM","Significant")]
    }else if(antibody == "ATAC"){
      active_mark <- read.csv(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_diff_in_20000_redundant_tads.csv"),row.names = 1)
      active_mark <- active_mark[,c("Geneid","LogFC.old.young","logCPM","Significant")]
    }else if(antibody == "RNA"){
      active_mark <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_in_20000_redundant_tads.csv"))
      active_mark <- active_mark[,c("X","logFC" ,"logCPM","Significant")]
      colnames(active_mark) <- c("Geneid","LogFC.old.young","logCPM","Significant")
    }
    colnames(active_mark)[4] <- "active_mark_significant"
    
    tad <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff_larger_250000.csv"))
    tad <- tad %>%
      separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)
    tad$start <- tad$start-1
    tad$Geneid <- paste0(tad$chr,":",tad$start,"-",tad$end)
    tad <- tad[,c("Geneid","Significant")]
    
    df <- merge(tad,active_mark,by="Geneid")
    to_plot <- df 
    to_plot$condition <- factor(to_plot$Significant, levels=c("Up","Stable","Down"))
    to_plot$tissue <- tissue_label_change(tissue)
    summary <- rbind(summary,to_plot)
  }
  to_plot <- summary %>%  
    group_by(tissue, condition) %>%  
    summarise(  
      median_logFC = median(LogFC.old.young, na.rm = TRUE),  
      count = n()  
    ) %>%   
    mutate(median_logFC = ifelse(count < 10, NA, median_logFC))
  
  p_value_summary <- data.frame(
    'Up' = rep(NA, length(tissues)), 
    'Stable' = rep(NA, length(tissues)), 
    'Down' = rep(NA, length(tissues))
  )
  rownames(p_value_summary) <- sapply(tissues, tissue_label_change)
  
  for(tissue in tissues){
    t_Up_summary <-  summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="Up"),]
    t_Stable_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="Stable"),]
    t_Down_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="Down"),]
    if(antibody %in% c("H3K27ac","H3K4me1","H3K4me3","RNA","ATAC","H3K36me3")){
      if(nrow(t_Up_summary) > 10){
        test <- wilcox.test(t_Up_summary$LogFC.old.young, t_Stable_summary$LogFC.old.young,alternative = "greater")
        p_value_summary[tissue_label_change(tissue),"Up"] <- test$p.value
      }
      if(nrow(t_Down_summary) > 10){
        test <- wilcox.test(t_Down_summary$LogFC.old.young,t_Stable_summary$LogFC.old.young,alternative = "less")
        p_value_summary[tissue_label_change(tissue),"Down"] <- test$p.value
      }
    }else{
      if(nrow(t_Up_summary) > 10){
        test <- wilcox.test(t_Up_summary$LogFC.old.young, t_Stable_summary$LogFC.old.young,alternative = "less")
        p_value_summary[tissue_label_change(tissue),"Up"] <- test$p.value
      }
      if(nrow(t_Down_summary) > 10){
        test <- wilcox.test(t_Down_summary$LogFC.old.young,t_Stable_summary$LogFC.old.young,alternative = "greater")
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
  
  p_value_summary$tissue <- rownames(p_value_summary)
  p_value_long <- reshape2::melt(p_value_summary,id.vars = "tissue")
  
  names(p_value_long) <- c("Tissue", "Type", "Label")
  colnames(to_plot)[c(1:3)] <- c("Tissue","Type","Value")
  merged_data <- merge(to_plot, p_value_long, by = c("Tissue", "Type"), all = TRUE)
  merged_data$Value[which(merged_data$Value > 1)] <- 1
  merged_data$Value[which(merged_data$Value < -1)] <- -1
  
  p <- ggplot(merged_data[which(merged_data$Type %in% c("Up","Down")),], aes(x = Type, y = Tissue, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         limits = c(-1, 1), midpoint = 0) +
    theme_minimal() +
    ggtitle(antibody)+
    geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
  ggsave(paste0("result/Sup_figures/TAD_",antibody,"_change_heatmap.pdf"),p,width = 6,height = 8)
}


