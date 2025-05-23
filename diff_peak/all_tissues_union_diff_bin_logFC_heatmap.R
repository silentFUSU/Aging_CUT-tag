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
antibody <- "H3K27me3"
common_increase <- read.csv(paste0("data/samples/all/",antibody,"/common_increase_10kb_bins_after_remove_batch_effect.csv"))
common_decrease <- read.csv(paste0("data/samples/all/",antibody,"/common_decrease_10kb_bins_after_remove_batch_effect.csv"))
# regions <- union(common_increase$Geneid[which(common_increase$n>8)],common_decrease$Geneid[which(common_decrease$n>8)])
# regions <- common_decrease$Geneid[which(common_decrease$n>20)]
# regions <- common_increase$Geneid[which(common_increase$n>20)]
regions <- c(common_decrease$Geneid[which(common_decrease$n>0)], common_increase$Geneid[which(common_increase$n>0)])
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()

for(tissue in tissues){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  df$LogFC.old.young[which(df$Significant=="Stable")] <- 0
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
k <- 5
kmeans_result <- kmeans(diff_summary, centers=k)
diff_summary$cluster <- kmeans_result$cluster  
diff_summary_sorted <- diff_summary[order(diff_summary$cluster),]  
data_for_heatmap <- diff_summary_sorted[, -ncol(diff_summary_sorted)]  

annotation <- diff_summary[,"cluster",drop=F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,breaks = breaks, annotation_row = annotation, color = color_palette)
# pheatmap::pheatmap(data_for_heatmap,cluster_rows = T,show_rownames = F,breaks = seq(-2, 2, length.out = 101), color = color_palette,)

df <- read.csv("data/samples/mammarygland/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv")
df <- df[which(df$Significant == "Down"),]
regions <- df$Geneid
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()
antibody <- "H3K9me3"
for(tissue in tissues){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
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
diff_summary <- diff_summary[order(diff_summary$`Mammary Gland`),]
pheatmap::pheatmap(diff_summary,cluster_rows = F,show_rownames = F,breaks = seq(-1, 1, length.out = 101) )


df <- read.csv("data/samples/ileum/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv")
df <- df[which(df$Significant=="Down"),]
regions <- df$Geneid
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()
antibody <- "H3K9me3"
for(tissue in tissues){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
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
diff_summary <- diff_summary[order(diff_summary$Ileum),]
pheatmap::pheatmap(diff_summary,cluster_rows = F,show_rownames = F,breaks = seq(-1, 1, length.out = 101) )
