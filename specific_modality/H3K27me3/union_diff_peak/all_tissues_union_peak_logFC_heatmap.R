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
antibody <- "H3K27me3"
tissue_num_summary <- read.csv("data/samples/all/H3K27me3/peaks_merged/merged_peaks_tissue_num.csv")
regions <- tissue_num_summary$label[which(tissue_num_summary$num >=27)]
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
diff_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_0_tissues_diff.csv"))
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
breaks <- c(seq(-2, -0.51, length.out = 40), seq(-0.5, 0.5, length.out = 20), seq(0.51, 2, length.out = 40))
p <- pheatmap::pheatmap(diff_summary,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")

diff_summary_rank <- diff_summary
diff_summary_rank$row_mean <- rowMeans(diff_summary_rank, na.rm = TRUE)
col_means <- colMeans(as.matrix(diff_summary_rank[, -ncol(diff_summary_rank)]), na.rm = TRUE)
diff_summary_rank_sorted_rows <- diff_summary_rank %>%
  arrange(row_mean)
diff_summary_rank_sorted_rows <- diff_summary_rank_sorted_rows[ , -ncol(diff_summary_rank_sorted_rows)]
diff_summary_rank_sorted <- diff_summary_rank_sorted_rows[, order(col_means)]
diff_summary_rank <- diff_summary_rank_sorted
p <- pheatmap::pheatmap(diff_summary_rank,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")



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

wcss <- sapply(1:20, function(k) {
  set.seed(1)
  kmeans_result <- kmeans(pca_df, centers = k, nstart = 25)
  return(kmeans_result$tot.withinss)
})
plot(1:20, wcss, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares")

silhouette_scores <- sapply(2:20, function(k) {
  cluster_assignment <- kmeans(pca_df, centers = k, nstart = 25)$cluster
  ss <- silhouette(cluster_assignment, dist(pca_df))
  return(mean(ss[, 3]))  
})
plot(2:20, silhouette_scores, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width")

set.seed(1)
umap_result <- umap::umap(as.matrix(pca_df))
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], names = rownames(pca_df))
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = names)) +
  geom_point() +
  theme_minimal() +
  labs(title = "UMAP Projection")

set.seed(1)
k <- 2
kmeans_df <- diff_summary
kmeans_df[is.na(kmeans_df)] <- 0
kmeans_result <- kmeans(kmeans_df, centers=k)
kmeans_result <- as.data.frame(kmeans_result$cluster)
colnames(kmeans_result) <- "cluster"
umap_df <- merge(umap_df,kmeans_result,by="row.names")
umap_df$cluster <- as.character(umap_df$cluster)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "UMAP Projection")

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
annotation_color <- list(cluster=setNames(c("#F8766D",  "#00BFC4"),c(1:2)))
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, color = color_palette, clustering_distance_cols="manhattan")
pheatmap::pheatmap(diff_summary,show_rownames = F,breaks = breaks, color = color_palette,  annotation_row = annotation,annotation_colors = annotation_color,clustering_distance_cols="manhattan",clustering_distance_rows="manhattan",clustering_method = "average")
dir.create("data/samples/all/H3K27me3/peaks_merged_exist_in_27_tissues_edger/")
write.csv(annotation,"data/samples/all/H3K27me3/peaks_merged_exist_in_27_tissues_edger/kmeans_annotation.csv")
dir.create("data/samples/all/H3K27me3/peaks_merged_exist_in_27_tissues_edger/bed/")
split_names <- strsplit(rownames(annotation), "[:-]")
annotation_df <- data.frame(
  chr = sapply(split_names, "[", 1),
  start = sapply(split_names, "[", 2),
  end = sapply(split_names, "[", 3),
  cluster = annotation$cluster
)
for(kmean in c(1:k)){
  write.table(annotation_df[which(annotation_df$cluster==kmean),c(1:3)],paste0("data/samples/all/H3K27me3/peaks_merged_exist_in_27_tissues_edger/bed/kmeans",kmean,"_uinon_recursion_peaks.bed"),col.names = F,row.names = F,append = F,quote = F,sep = "\t")
}
