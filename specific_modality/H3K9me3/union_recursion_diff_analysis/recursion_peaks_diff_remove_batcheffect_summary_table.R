rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(stringr)
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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
antibody <- "H3K9me3"
window_size="5000"
gap_size="10000"
summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion_diff_after_remove_batch_effect.csv"))
  t_summary <- data.frame(tissue=tissue_label_change(tissue),
                          all_peaks=nrow(df),
                          increased=nrow(df[which(df$Significant=="Up"),]),
                          decreased=nrow(df[which(df$Significant=="Down"),]),
                          long_peaks=nrow(df[which(df$Length >= 200000),]),
                          increased_long=nrow(df[which(df$Length >= 200000 & df$Significant=="Up"),]),
                          decreased_long=nrow(df[which(df$Length >= 200000 & df$Significant=="Down"),]))  
  summary <- rbind(summary,t_summary)
  }

write.csv(summary,"data/samples/all/H3K9me3/recursion_peaks_diff_table/all_tissue_diff_recursion_peaks_summary.csv")


