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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()
diff_summary_Discrete <- data.frame()
antibody <- "H3K9me3"
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Length >=200000),c("Geneid","LogFC.old.young","Significant")]
  df_Discrete <- df[,-2]
  df_Discrete$value <- 0
  df_Discrete$value[which(df_Discrete$Significant=="Down")] <- -1
  df_Discrete$value[which(df_Discrete$Significant=="Up")] <- 1
  df_Discrete <- df_Discrete[,-2]
  colnames(df_Discrete)[2] <- tissue_label_change(tissue)
  if(nrow(diff_summary_Discrete) == 0){
    diff_summary_Discrete <- df_Discrete
  }else{
    diff_summary_Discrete <- merge(diff_summary_Discrete,df_Discrete,by="Geneid",all=T)    
  }
  
  df <- df[,-3]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(diff_summary) == 0){
    diff_summary <- df
  }else{
    diff_summary <- merge(diff_summary,df,by="Geneid",all=T)    
  }
}
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
colnames(annotation)[1] <- "Geneid"

diff_summary <- merge(diff_summary,annotation,by="Geneid",all=T)
rownames(diff_summary) <- diff_summary$Geneid
diff_summary$cluster[is.na(diff_summary$cluster)] <- "Stable"
to_plot <- diff_summary[order(diff_summary$cluster),]
to_plot <- to_plot[,-c(1,ncol(to_plot))]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
annotation_color <- list(cluster=setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(1:4,"Stable")))
annotation <- data.frame(Geneid = diff_summary$Geneid,cluster=diff_summary$cluster)
rownames(annotation) <- annotation$Geneid
annotation <- annotation[,-1,drop=F]
tissue_order <- c("Kidney","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","IWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows = F,show_rownames = F,
                   breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, 
                   color = color_palette, clustering_distance_cols="manhattan",
                   border_color = "black",gaps_row = c(283,283+271,283+271+347,283+271+347+187),cluster_cols = F)

diff_summary <- diff_summary[,-1]
write.csv(diff_summary[order(diff_summary$cluster),c(tissue_order,"cluster")],"data/samples/all/H3K9me3/recursion_peaks_diff_table/all_tissue_diff_recursion_peaks_heatmap_matrix.csv")

#####
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
colnames(annotation)[1] <- "Geneid"
diff_summary_Discrete <- merge(diff_summary_Discrete,annotation,by="Geneid",all=T)
rownames(diff_summary_Discrete) <- diff_summary_Discrete$Geneid
diff_summary_Discrete$cluster[is.na(diff_summary_Discrete$cluster)] <- "Stable"
to_plot <- diff_summary_Discrete[order(diff_summary_Discrete$cluster),]
to_plot <- to_plot[,-c(1,ncol(to_plot))]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
annotation_color <- list(cluster=setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(1:4,"Stable")))
annotation <- data.frame(Geneid = diff_summary_Discrete$Geneid,cluster=diff_summary_Discrete$cluster)
rownames(annotation) <- annotation$Geneid
annotation <- annotation[,-1,drop=F]
tissue_order <- c("Kidney","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","IWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows = F,show_rownames = F,
                   breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, 
                   color = color_palette, clustering_distance_cols="manhattan",
                   border_color = "black",gaps_row = c(283,283+271,283+271+347,283+271+347+187),cluster_cols = F)

diff_summary <- diff_summary[,-1]
