rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
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

tissue <- "lung"
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_num <-27
rpkm_log2FC_summary <- data.frame()
condition <- "peak"
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(condition=="domain"){
    tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts.summary"),header = T,row.names = 1)
  }else{
    tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.counts.summary"),header = T,row.names = 1)
  }

  if(tissue %in% c("mammarygland","uterus","ovary")){
    tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
  }
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  colnames(summary) <- gsub(pattern,"\\1",colnames(summary))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  summary <- summary[-2,search_table$sample_name]
  total_reads <- colSums(summary)
  length <- as.numeric(tab$Length)
  
  rpkm <- sweep(counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  
  rpkm_young <- rpkm[,search_table$sample_name[which(search_table$age=="3m")]]
  rpkm_old <- rpkm[,search_table$sample_name[which(search_table$age=="24m")]]
  
  rpkm_young$mean_young <- rowMeans(rpkm_young)
  rpkm_old$mean_old <- rowMeans(rpkm_old)
  
  rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
  rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
  colnames(rpkm_mean_summary)[which(colnames(rpkm_mean_summary)=="log2FC")] <- tissue_label_change(tissue)
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  if(nrow(rpkm_log2FC_summary )==0){
    rpkm_log2FC_summary <- rpkm_mean_summary[,c(1,4)]
  }else{
    rpkm_log2FC_summary <- merge(rpkm_log2FC_summary,rpkm_mean_summary[,c(1,4)],by="Geneid",all=T)
  }
}

rownames(rpkm_log2FC_summary) <- rpkm_log2FC_summary$Geneid
rpkm_log2FC_summary <- rpkm_log2FC_summary[,-1]
rpkm_log2FC_summary$row_mean <- rowMeans(rpkm_log2FC_summary, na.rm = TRUE)
col_means <- colMeans(as.matrix(rpkm_log2FC_summary[, -ncol(rpkm_log2FC_summary)]), na.rm = TRUE)
rpkm_log2FC_summary_sorted_rows <- rpkm_log2FC_summary %>%
  arrange(row_mean)
rpkm_log2FC_summary_sorted_rows <- rpkm_log2FC_summary_sorted_rows[ , -ncol(rpkm_log2FC_summary_sorted_rows)]
rpkm_log2FC_summary_sorted <- rpkm_log2FC_summary_sorted_rows[, order(col_means)]
rpkm_log2FC_summary <- rpkm_log2FC_summary_sorted

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-5, -0.51, length.out = 40), seq(-0.5, 0.5, length.out = 20), seq(0.51, 5, length.out = 40))
p <- pheatmap::pheatmap(rpkm_log2FC_summary,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")
p <- pheatmap::pheatmap(rpkm_log2FC_summary,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")




pca_df <- rpkm_log2FC_summary
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
kmeans_df <- rpkm_log2FC_summary
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


to_plot <- merge(rpkm_log2FC_summary,kmeans_result,by="row.names") 
rpkm_log2FC_summary_sorted <- to_plot[order(to_plot$cluster),]
rownames(rpkm_log2FC_summary_sorted) <- rpkm_log2FC_summary_sorted$Row.names
rpkm_log2FC_summary_sorted <- rpkm_log2FC_summary_sorted[-which(colnames(rpkm_log2FC_summary_sorted) %in% c("Row.names","cluster"))]
data_for_heatmap <- rpkm_log2FC_summary_sorted
annotation <- to_plot[,c("Row.names","cluster"),drop=F]
rownames(annotation) <- annotation$Row.names
annotation <- annotation[,-1,drop=F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-5, -0.51, length.out = 40), seq(-0.5, 0.5, length.out = 20), seq(0.51, 5, length.out = 40))
annotation_color <- list(cluster=setNames(c("#F8766D", "#7CAE00", "#00BFC4")[1:k],c(1:2)))
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,show_rownames = F,breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color, color = color_palette, clustering_distance_cols="manhattan")
pheatmap::pheatmap(rpkm_log2FC_summary,show_rownames = F,breaks = breaks, color = color_palette,  annotation_row = annotation,annotation_colors = annotation_color,clustering_distance_cols="manhattan",clustering_distance_rows="manhattan",clustering_method = "average")
if(condition=="domain"){
  dir.create("data/samples/all/H3K27me3/edd_domain_merged/")
  write.csv(annotation,"data/samples/all/H3K27me3/edd_domain_merged/kmeans_annotation_RPKM.csv")
}else{
  dir.create("data/samples/all/H3K27me3/peaks_merged/")
  write.csv(annotation,paste0("data/samples/all/H3K27me3/peaks_merged/kmeans_annotation_larger_",tissue_num,"_tissues_RPKM.csv"))
}

split_names <- strsplit(rownames(annotation), "[:-]")
annotation_df <- data.frame(
  chr = sapply(split_names, "[", 1),
  start = sapply(split_names, "[", 2),
  end = sapply(split_names, "[", 3),
  cluster = annotation$cluster
)
for(kmean in c(1:k)){
  if(condition == "domain"){
    write.table(annotation_df[which(annotation_df$cluster==kmean),c(1:3)],paste0("data/samples/all/H3K27me3/edd_domain_merged/kmeans",kmean,"_edd_domain_merged_RPKM.bed"),col.names = F,row.names = F,append = F,quote = F,sep = "\t")
  }else{
    write.table(annotation_df[which(annotation_df$cluster==kmean),c(1:3)],paste0("data/samples/all/H3K27me3/peaks_merged/kmeans",kmean,"_larger_",tissue_num,"_tissues_RPKM.bed"),col.names = F,row.names = F,append = F,quote = F,sep = "\t")
  }
}
