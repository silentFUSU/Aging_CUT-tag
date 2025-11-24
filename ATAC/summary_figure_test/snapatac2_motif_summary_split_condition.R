rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
extract_before_bracket <- function(s) {  
  parts <- strsplit(s, "\\(")[[1]]  
  return(parts[1])  
} 
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

motif_data_frame <- list(up=data.frame(),down=data.frame())
conditions <- c("up","down")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
data_path <- "data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg/"
if(data_path=="data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg/"){
  label <- "(using stable peaks as background)"
}
increase_count <- read.csv("data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg/increase_common_motif_count.csv")
decrease_count <- read.csv("data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg/decrease_common_motif_count.csv")
common_count <- list(up=increase_count,down=decrease_count)
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(conditions))){
    condition <- conditions[j]
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_",str_to_title(condition),"_sorted.bed.csv"))){
      motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_",str_to_title(condition),"_sorted.bed.csv"))
      motif$adjusted.p.value[which(motif$log2.fold.change. < 0 )] <- 1
      motif$id <- ifelse(
        grepl("\\(.*\\)", motif$id),
        sub(".*\\((.*?)\\).*", "\\1", motif$id), 
        sub("^M\\d+_2\\.00\\s*", "", motif$id) 
      )
      
      colnames(motif)[which(colnames(motif)=="adjusted.p.value")] <- tissue_label_change(tissue)
      if(i == 1){
        motif_data_frame[[condition]] <- motif[,c("id",tissue_label_change(tissue))]
      }else{
        motif_data_frame[[condition]] <- merge(motif_data_frame[[condition]],motif[,c("id",tissue_label_change(tissue))],by="id",all=T)
        }
      }
    }
  }

condition <- "down"
to_plot <- motif_data_frame[[condition]]
rownames(to_plot) <- to_plot$id
to_plot <- to_plot[,-1]
replace_zeros <- function(column) {
  non_zero_min <- min(column[column != 0], na.rm = TRUE)
  if (is.infinite(non_zero_min)) {
    non_zero_max <- 0
  }
  column[column == 0] <- non_zero_min
  return(column)
}
to_plot <- as.data.frame(apply(to_plot, 2, replace_zeros))
to_plot <- as.data.frame(apply(to_plot, 2, function(column) -log10(column)))
# columns_to_keep <- apply(to_plot, 2, function(column) any(column >= -log10(0.05)))
# to_plot <- to_plot[, columns_to_keep]

if(condition =="up"){
  color <- c("#ffe2e2","red")
  clusters <- c("cluster1","cluster2","tissue specific")
  threshold <- 5
}else{
  color <- c("#defcf9","blue")
  clusters <- c("larger6","cluster2","tissue specific")
  threshold <- 6
}
color_palette <- c(
  rep("white", 20), 
  colorRampPalette(color)(80) 
)
breaks <- c(seq(0, -log10(0.05), length.out = 20), seq(-log10(0.05)+0.0001, 3, length.out = 80))
to_plot <- as.data.frame(to_plot[common_count[[condition]]$id,])
p <- pheatmap::pheatmap(to_plot, cluster_rows =T,cluster_cols = T,color = color_palette,breaks=breaks,show_rownames = F,clustering_distance_cols = "euclidean")
tissue_order <- p$gtable$grobs[[4]]$label
to_plot <- to_plot[,tissue_order]

common_count_annotation <- common_count[[condition]]
common_count_annotation$categories <- clusters[1]
common_count_annotation$categories[which(common_count_annotation$n >=2 & common_count_annotation$n < threshold)] <- clusters[2]
common_count_annotation$categories[which(common_count_annotation$n ==1)] <- clusters[3]
common_count_annotation$categories <- factor(common_count_annotation$categories, levels=clusters)
common_count_annotation <- common_count_annotation[order(common_count_annotation$categories),]
to_plot <- to_plot[common_count_annotation$id,]
gaps_row <- cumsum(table(common_count_annotation$categories))
annotation_row <-common_count_annotation[,c("id","categories")]
rownames(annotation_row) <- annotation_row$id
annotation_row <- annotation_row[,-1,drop=F]

pheatmap::pheatmap(
  to_plot, 
  annotation_row = annotation_row, 
  cluster_rows = F,  
  clustering_method = "complete", 
  gaps_row = gaps_row,
  breaks = breaks,
  color = color_palette,
  border_color = "black",show_rownames = F
)
cluster1 <- common_count[[condition]]$id[which(common_count[[condition]]$n >=threshold)]
to_plot_cluster1 <- to_plot[cluster1,]
p1 <- pheatmap::pheatmap(to_plot_cluster1, cluster_rows =T,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 7)
p1_label <- p1$gtable$grobs[[4]]$label
rowmeans <- as.data.frame(rowMeans(to_plot_cluster1))
colnames(rowmeans) <- "means"
rowmeans <- rowmeans[order(rowmeans$means,decreasing = T),,drop=F]
to_plot_cluster1 <- to_plot_cluster1[rownames(rowmeans),]

colmeans <- as.data.frame(colMeans(to_plot_cluster1))
colnames(colmeans) <- "means"
colmeans  <- colmeans [order(colmeans$means,decreasing = T),,drop=F]
to_plot_cluster1 <- to_plot_cluster1[,rownames(colmeans)]
pheatmap::pheatmap(to_plot_cluster1, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 7)

# columns_to_keep <- apply(to_plot, 2, function(column) sum(column >= -log10(0.05)) >= 10)
# to_plot_cluster1 <- to_plot_cluster1[,columns_to_keep]
# pheatmap::pheatmap(to_plot_cluster1, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 7)


cluster2 <- common_count[[condition]]$id[which(common_count[[condition]]$n >=2 & common_count[[condition]]$n < threshold)]
to_plot_cluster2 <- to_plot[cluster2,]
p2 <- pheatmap::pheatmap(to_plot_cluster2, cluster_rows =T,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T)
p2_label <- p2$gtable$grobs[[4]]$label

tissue_specific <- common_count[[condition]]$id[which(common_count[[condition]]$n ==1)]
to_plot_tissue_specific <- to_plot[tissue_specific,]
p3 <- pheatmap::pheatmap(to_plot_tissue_specific, cluster_rows =T,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T)
p3_label <- p3$gtable$grobs[[4]]$label

motif_order <- c(p1_label,p2_label,p3_label)
pheatmap::pheatmap(
  to_plot[motif_order,], 
  annotation_row = annotation_row, 
  cluster_rows = F,  
  clustering_method = "complete", 
  gaps_row = gaps_row,
  breaks = breaks,
  color = color_palette,
  border_color = "black",
  show_rownames = F
)
