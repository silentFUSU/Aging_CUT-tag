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

regions <- unique(c(common_decrease$Geneid[which(common_decrease$n>0)], common_increase$Geneid[which(common_increase$n>0)]))
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
    diff_summary <- merge(diff_summary,df,by="Geneid",all=T)    
  }
}

rownames(diff_summary) <- diff_summary$Geneid
diff_summary <- diff_summary[,-1]

pca_df <- diff_summary
pca_df[is.na(pca_df)] <- 0
pca <- prcomp(pca_df)
to_plot <- data.frame(pca$x)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2)) + 
  geom_point(size=1) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))


### test kmeans count
wcss <- sapply(1:20, function(k) {
  set.seed(1)
  kmeans_result <- kmeans(pca_df, centers = k, nstart = 25)
  return(kmeans_result$tot.withinss)
})
plot(1:20, wcss, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares")

### kmeans clustering
set.seed(1)
k <- 4
kmeans_df <- diff_summary
kmeans_df[is.na(kmeans_df)] <- 0
kmeans_result <- kmeans(kmeans_df, centers=k)
kmeans_result <- as.data.frame(kmeans_result$cluster)
colnames(kmeans_result) <- "cluster"

to_plot <- merge(diff_summary,kmeans_result,by="row.names") 
diff_summary_sorted <- to_plot[order(to_plot$cluster),]
rownames(diff_summary_sorted) <- diff_summary_sorted$Row.names
diff_summary_sorted <- diff_summary_sorted[-which(colnames(diff_summary_sorted) %in% c("Row.names","cluster"))]
data_for_heatmap <- diff_summary_sorted
annotation <- to_plot[,c("Row.names","cluster"),drop=F]
rownames(annotation) <- annotation$Row.names
annotation <- annotation[,-1,drop=F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
annotation_color <- list(cluster=setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),c(1:4)))
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, color = color_palette, clustering_distance_cols="manhattan")
write.csv(annotation,"data/samples/all/H3K27me3/recursion_bin_diff_table/kmeans_annotation.csv")
