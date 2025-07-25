rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(factoextra)
library(cluster)
library(umap)
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

antibody <- "H3K9me3"
common_decrease <- read.csv(paste0("data/samples/all/",antibody,"/common_decrease-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_decrease <- common_decrease[which((common_decrease$end - common_decrease$start +1) > 200000),]
regions <- unique(c(common_decrease$Geneid[which(common_decrease$n>0)]))
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Geneid %in% regions),c("Geneid","LogFC.old.young","Significant")]
  df <- df[,-3]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(diff_summary) == 0){
    diff_summary <- df
  }else{
    diff_summary <- merge(diff_summary,df,by="Geneid",all=T)    
  }
}

rownames(diff_summary) <- diff_summary$Geneid
diff_summary <- diff_summary[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
rownames(annotation) <- annotation$X
annotation$cluster <- as.character(annotation$cluster)
annotation <- annotation[,-1,drop=F]
p <- pheatmap::pheatmap(diff_summary,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan",annotation_row = annotation)


