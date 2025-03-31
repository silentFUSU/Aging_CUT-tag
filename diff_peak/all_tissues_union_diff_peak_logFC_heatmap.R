rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
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
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
# common_increase <- read.csv("data/samples/all/H3K9me3/common_increase_union_peaks_after_remove_batch_effect.csv")
# common_decrease <- read.csv("data/samples/all/H3K9me3/common_decrease_union_peaks_after_remove_batch_effect.csv")
# regions <- union(common_increase$Geneid[which(common_increase$n>8)],common_decrease$Geneid[which(common_decrease$n>8)])
antibody <- "H3K27me3"
common_increase <- read.csv(paste0("data/samples/all/",antibody,"/common_increase-W5000-G10000-E100_union_peaks_after_remove_batch_effect.csv"))
common_increase <- common_increase[which((common_increase$end - common_increase$start +1) > 100000),]
common_decrease <- read.csv(paste0("data/samples/all/",antibody,"/common_decrease-W5000-G10000-E100_union_peaks_after_remove_batch_effect.csv"))
common_decrease <- common_decrease[which((common_decrease$end - common_decrease$start +1) > 100000),]
regions <- common_increase$Geneid[which(common_increase$n>5)]
# regions <- common_increase$Geneid[which(common_increase$n>5)]
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()

for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W5000-G10000-E100_bedtools_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Geneid %in% regions),c("Geneid","LogFC.old.young")]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(diff_summary) == 0){
    diff_summary <- df
  }else{
    diff_summary <- merge(diff_summary,df,by="Geneid")    
  }
}
rownames(diff_summary) <- diff_summary$Geneid
diff_summary <- diff_summary[,-1]
set.seed(1)
k <- 10
kmeans_result <- kmeans(diff_summary, centers=k)  
diff_summary$cluster <- kmeans_result$cluster  
diff_summary_sorted <- diff_summary[order(diff_summary$cluster),]  
data_for_heatmap <- diff_summary_sorted[, -ncol(diff_summary_sorted)]  
distance <- dist(data_for_heatmap)
hc_rows <- hclust(distance)
pheatmap::pheatmap(data_for_heatmap,cluster_rows = T,show_rownames = F,breaks = seq(-1, 1, length.out = 101))

common_increase$Length <- common_increase$end - common_increase$start + 1
common_decrease$Length <- common_decrease$end - common_decrease$start + 1
common_increase <- common_increase[which(common_increase$n > 10),]
common_decrease <- common_decrease[which(common_decrease$n > 10),]
common_decrease$condition <- "Common decrease"
common_increase$condition <- "Common increase"
length_summary <- rbind(common_decrease[,c("Length","condition")],common_increase[,c("Length","condition")])
t.test(common_increase$Length,common_decrease$Length)
