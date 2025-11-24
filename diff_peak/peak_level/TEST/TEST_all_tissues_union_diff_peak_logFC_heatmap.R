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
antibody <- "H3K9me3"
common_increase <- read.csv(paste0("data/samples/all/",antibody,"/common_increase-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_increase <- common_increase[which((common_increase$end - common_increase$start +1) > 200000),]
common_decrease <- read.csv(paste0("data/samples/all/",antibody,"/common_decrease-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_decrease <- common_decrease[which((common_decrease$end - common_decrease$start +1) > 200000),]
# regions <- common_increase$Geneid[which(common_increase$n>0)]
# regions <- common_decrease$Geneid[which(common_decrease$n>0)]
regions <- c(common_increase$Geneid[which(common_increase$n>0)],common_decrease$Geneid[which(common_decrease$n>0)])
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()

for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Geneid %in% regions),c("Geneid","LogFC.old.young","Significant")]
  # df$LogFC.old.young[which(df$Significant=="Stable")] <- 0
  df <- df[,-3]
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
k <- 4
kmeans_result <- kmeans(diff_summary, centers=k)
diff_summary$cluster <- kmeans_result$cluster  
diff_summary_sorted <- diff_summary[order(diff_summary$cluster),]  
data_for_heatmap <- diff_summary_sorted[, -ncol(diff_summary_sorted)]  


annotation <- diff_summary[,"cluster",drop=F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,breaks = breaks, annotation_row = annotation, color = color_palette, clustering_distance_cols="manhattan")
write.csv(annotation,"data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
split_names <- strsplit(rownames(annotation), "[:-]")
annotation_df <- data.frame(
  chr = sapply(split_names, "[", 1),
  start = sapply(split_names, "[", 2),
  end = sapply(split_names, "[", 3),
  cluster = annotation$cluster
)
for(kmean in c(1:k)){
  write.table(annotation_df[which(annotation_df$cluster==kmean),c(1:3)],paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/kmeans",kmean,"_uinon_recursion_peaks.bed"),col.names = F,row.names = F,append = F,quote = F,sep = "\t")
}
df <-read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv")
df <- df[,-1]
df <- merge(df,annotation,by.x="Geneid",by.y="row.names",all=T)
write.csv(df,"data/samples/all/H3K9me3/recursion_peaks_diff_table/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv",row.names = F)

common_increase$Length <- common_increase$end - common_increase$start + 1
common_decrease$Length <- common_decrease$end - common_decrease$start + 1
common_increase <- common_increase[which(common_increase$n > 10),]
common_decrease <- common_decrease[which(common_decrease$n > 10),]
common_decrease$condition <- "Common decrease"
common_increase$condition <- "Common increase"
length_summary <- rbind(common_decrease[,c("Length","condition")],common_increase[,c("Length","condition")])
t.test(common_increase$Length,common_decrease$Length)
