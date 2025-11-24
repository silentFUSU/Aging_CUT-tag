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
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","MEF")
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
tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","IWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows = F,show_rownames = F,
                   breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, 
                   color = color_palette, clustering_distance_cols="manhattan",
                   border_color = "black",gaps_row = c(283,283+271,283+271+347,283+271+347+187),cluster_cols = F)


# annotation[rownames(annotation) %in% grep("chrY", rownames(annotation), value = TRUE) & annotation$cluster == 1, "cluster"] <- "kmeans1-chrY"

# kmeans1_chrY <- to_plot[rownames(annotation[which(annotation$cluster=="kmeans1-chrY"),,drop=F]),]
# kmeans1_chrY_medians <- apply(kmeans1_chrY, 2, median, na.rm = TRUE)
# kmeans1_chrY_medians <- data.frame(tissues=colnames(to_plot),kmeans1_chrY_median=kmeans1_chrY_medians)
# kmeans1_chrY_medians$kmeans1_chrY_median[which(kmeans1_chrY_medians$kmeans1_chrY_median >0)] <- 0

kmeans1 <- to_plot[rownames(annotation[which(annotation$cluster=="1"),,drop=F]),]
kmeans1_medians <- apply(kmeans1, 2, median, na.rm = TRUE)
kmeans1_medians <- data.frame(tissues=colnames(to_plot),kmeans1_median=kmeans1_medians)
kmeans1_medians$kmeans1_median[which(kmeans1_medians$kmeans1_median >0)] <- 0


kmeans3 <- to_plot[rownames(annotation[which(annotation$cluster=="3"),,drop=F]),]
kmeans3_medians <- apply(kmeans3, 2, median, na.rm = TRUE)
kmeans3_medians <- data.frame(tissues=colnames(to_plot),kmeans3_median=kmeans3_medians)
kmeans3_medians$kmeans3_median[which(kmeans3_medians$kmeans3_median >0 )] <- 0

tissue_order <- merge(kmeans1_medians,kmeans3_medians,by="tissues")
tissue_order$weighted_value <- tissue_order$kmeans1_median + -1* tissue_order$kmeans3_median 
tissue_order <- tissue_order[order(tissue_order$weighted_value),]
 ggplot(tissue_order, aes(x = reorder(tissues, weighted_value), y = weighted_value)) +
  geom_bar(stat = "identity") +
  labs(x = "Tissues", y = "Weighted Value", title = "Bar Plot of Tissues Weighted Value") +
  theme_bw()+
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 90, hjust = 1,vjust = 0.5),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 

p <-ggplot(tissue_order, aes(x = reorder(tissues, weighted_value), y = kmeans1_median)) +
  geom_bar(stat = "identity") +
  labs(x = "Tissues", y = "Weighted Value", title = "Bar Plot of Tissues Weighted Value") +
  theme_bw()+
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 90, hjust = 1,vjust = 0.5),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  )+ylim(-0.4,0)
ggsave("result/Sup_figures/kmeans1_logFC_bar.pdf",p,width = 8,height = 4)

p <- ggplot(tissue_order, aes(x = reorder(tissues, weighted_value), y = kmeans3_median)) +
  geom_bar(stat = "identity") +
  labs(x = "Tissues", y = "Weighted Value", title = "Bar Plot of Tissues Weighted Value") +
  theme_bw()+
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 90, hjust = 1,vjust = 0.5),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) +ylim(-0.4,0)
ggsave("result/Sup_figures/kmeans3_logFC_bar.pdf",p,width = 8,height = 4)

# kmeans1_chrY_medians$tissues <- factor(kmeans1_chrY_medians$tissues,levels = tissue_order$tissues)
# ggplot(kmeans1_chrY_medians, aes(x = tissues, y = kmeans1_chrY_median)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Tissues", y = "Weighted Value", title = "Bar Plot of Tissues Weighted Value") +
#   theme_bw()+
#   theme(  
#     axis.title.x = element_text(size = 14),  
#     axis.title.y = element_text(size = 14),  
#     axis.text.x = element_text(size = 14,angle = 90, hjust = 1,vjust = 0.5),  
#     axis.text.y = element_text(size = 14),  
#     plot.title = element_text(size = 16, face = "bold")
#   ) 


tissue_order <- tissue_order$tissues
# annotation$cluster[which(annotation$cluster=="kmeans1-chrY")] <- 1
to_plot <- to_plot[,tissue_order]
breaks <- c(seq(-1, -0.41, length.out = 40), seq(-0.4, 0.4, length.out = 20), seq(0.41, 1, length.out = 40))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap::pheatmap(to_plot,cluster_rows = F,show_rownames = F,
                   breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, 
                   color = color_palette,
                   border_color = "black",gaps_row = c(283,283+271,283+271+347,283+271+347+187),cluster_cols = F,filename = paste0("result/figures/H3K9me3_heatmap.pdf"),
                   width = 6,height =9)
