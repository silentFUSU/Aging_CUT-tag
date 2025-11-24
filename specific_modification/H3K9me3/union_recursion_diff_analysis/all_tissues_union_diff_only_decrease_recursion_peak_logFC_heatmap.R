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
umap_df$chromosome <- sub(":(.*)", "", umap_df$names)
umap_df$chromosome[which(umap_df$chromosome != "chrY")] <- "other"
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = names,color=chromosome)) +
  geom_point() +
  theme_minimal() +
  labs(title = "UMAP Projection")

### kmeans clustering
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
# color_palette <- colorRampPalette(c("#6a65d8","white","#e23e57"))(100) 
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
# annotation_color <- list(cluster=setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),c(1:4)))
annotation_color <- list(cluster=setNames(c("#f6416c", "#f8f3d4"),c(1:2)))
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,
                   breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, 
                   color = color_palette, clustering_distance_cols="manhattan",
                   border_color = "black",
                   width = 6,height =9) 

### tissues clustering
pca_df <- as.data.frame(t(diff_summary))
pca_df[is.na(pca_df)] <- 0
pca <- prcomp(pca_df)
to_plot <- data.frame(pca$x)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
ggplot(to_plot, aes(x=PC1, y=PC2)) + 
  geom_point(size=1) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))

wcss <- sapply(1:10, function(k) {
  set.seed(1)
  kmeans_result <- kmeans(pca_df, centers = k, nstart = 25)
  return(kmeans_result$tot.withinss)
})
plot(1:10, wcss, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares")

silhouette_scores <- sapply(2:10, function(k) {
  cluster_assignment <- kmeans(pca_df, centers = k, nstart = 25)$cluster
  ss <- silhouette(cluster_assignment, dist(pca_df))
  return(mean(ss[, 3]))  
})
plot(2:10, silhouette_scores, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width")

set.seed(1)
umap_result <- umap::umap(as.matrix(pca_df))
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], names = rownames(pca_df))
tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = names,color=names)) +
  geom_point(size=3) +
  scale_color_manual(values = color) +
  theme_minimal() +
  labs(title = "UMAP Projection")

set.seed(1)
k <- 3
kmeans_df_tissue <- as.data.frame(t(diff_summary))
kmeans_df_tissue[is.na(kmeans_df_tissue)] <- 0
kmeans_result_tissue <- kmeans(kmeans_df_tissue, centers=k)
kmeans_result_tissue <- as.data.frame(kmeans_result_tissue$cluster)
colnames(kmeans_result_tissue) <- "cluster"
umap_df <- merge(umap_df,kmeans_result_tissue,by="row.names")
umap_df$cluster <- as.character(umap_df$cluster)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title = "UMAP Projection")

kmeans_result_tissue <- kmeans_result_tissue[order(kmeans_result_tissue$cluster),,drop=F]
colnames(kmeans_result_tissue) <- "tissue"
kmeans_result_tissue$tissue <- as.character(kmeans_result_tissue$tissue)
data_for_heatmap <- data_for_heatmap[,rownames(kmeans_result_tissue)]
breaks <- c(seq(-1, -0.55, length.out = 40), seq(-0.5,0.5, length.out = 20), seq(0.55, 1, length.out = 40))
annotation_color <- list(cluster=setNames(c("#f6416c", "#f8f3d4"),c(1:2)), tissue=setNames(c("#ff165d","#00b8a9","#07689f"),c(1:3)))
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,cluster_cols = F,
                   breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, 
                   color = color_palette,annotation_col = kmeans_result_tissue,
                   border_color = "black")
