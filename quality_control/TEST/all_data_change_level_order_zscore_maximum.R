rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
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
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
### Histone modification
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
summary <- data.frame()
for(antibody in antibodys){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  summary_per_antibody <- data.frame()
  for(tissue in tissues){
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
    df <- df[which(df$Significant != "Stable"),]
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_in_young_old_merge-W1000-G3000-E100.bed"))    
    }else{
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_in_young_old_merge_macs_narrowpeak.bed"))
    }
    df <- df[which(df$Geneid %in% peaks$V4),]
    t_summary_per_antibody <- data.frame(tissue=tissue_label_change(tissue),count=nrow(df))
    summary_per_antibody <- rbind(summary_per_antibody,t_summary_per_antibody)
  }
  colnames(summary_per_antibody)[2] <- antibody
  summary_per_antibody <- summary_per_antibody[,c(1,2)]
  if(nrow(summary)==0){
    summary <- summary_per_antibody
  }else{
    summary <- merge(summary,summary_per_antibody,by="tissue")
  }
}

##### gene_expression
summary_RNA <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  df <- df[which(df$Significant !="Stable"),]
  t_summary_RNA <- data.frame(tissue=tissue_label_change(tissue),count=nrow(df))  
  summary_RNA <- rbind(summary_RNA,t_summary_RNA)
}
colnames(summary_RNA)[2] <- "RNA"
summary <- merge(summary,summary_RNA,by="tissue")

##### DNA methylation
summary_DNA <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01.txt"))
  t_summary_DNA <- data.frame(tissue=tissue_label_change(tissue),count=nrow(df))
  summary_DNA <- rbind(summary_DNA,t_summary_DNA)
}
colnames(summary_DNA)[2] <- "WGBS"
summary <- merge(summary,summary_DNA,by="tissue")

##### ATAC 
summary_ATAC <- data.frame()
for(tissue in tissues){
  df_up <- read.csv(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/up/",tissue,"_Up.bed"))
  df_down <- read.csv(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/down/",tissue,"_Down.bed"))
  t_summary_ATAC <- data.frame(tissue=tissue_label_change(tissue),count=(nrow(df_down)+nrow(df_up)))
  summary_ATAC <- rbind(summary_ATAC,t_summary_ATAC)
}
colnames(summary_ATAC)[2] <- "ATAC"
summary <- merge(summary,summary_ATAC,by="tissue")
rownames(summary) <- summary$tissue

#### zscore
z_scores <- as.data.frame(lapply(summary[,-1], scale))
df_z_scores <- cbind(tissue = summary$tissue, z_scores)
rownames(df_z_scores) <- df_z_scores$tissue
to_plot <- df_z_scores[,-1]
to_plot_rank <- to_plot
to_plot_rank$means <- rowMeans(to_plot_rank)
to_plot_rank <- to_plot_rank[order(to_plot_rank$means,decreasing = T),]
to_plot <- to_plot[rownames(to_plot_rank),]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-4, -1.1, length.out = 40), seq(-1, 1, length.out = 20), seq(1.1, 4, length.out = 40))
pheatmap::pheatmap(to_plot,cluster_rows = F,color = color_palette,breaks = breaks,cluster_cols = F)

summary_male <- summary[-which(summary$tissue %in% c("Mammary Gland","Uterus","Ovary")),]
z_scores <- as.data.frame(lapply(summary_male[,-1], scale))
df_z_scores <- cbind(tissue = summary_male$tissue, z_scores)
rownames(df_z_scores) <- df_z_scores$tissue
to_plot <- df_z_scores[,-1]
to_plot_rank <- to_plot
to_plot_rank$means <- rowMeans(to_plot_rank)
to_plot_rank <- to_plot_rank[order(to_plot_rank$means,decreasing = T),]
to_plot <- to_plot[rownames(to_plot_rank),]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-4, -1.1, length.out = 40), seq(-1, 1, length.out = 20), seq(1.1, 4, length.out = 40))
pheatmap::pheatmap(to_plot,cluster_rows = F,color = color_palette,breaks = breaks,cluster_cols = F)

#### maximum
to_plot <- summary
to_plot <- to_plot[-which(to_plot$tissue %in% c("Mammary Gland","Uterus","Ovary")),]
to_plot <- to_plot[,-1]
max_values <- apply(to_plot, 2, max)
normalized_to_plot <- sweep(to_plot, 2, max_values, FUN = "/")
to_plot_rank <- to_plot
to_plot_rank$means <- rowMeans(to_plot_rank)
to_plot_rank <- to_plot_rank[order(to_plot_rank$means,decreasing = T),]
normalized_to_plot <- normalized_to_plot[rownames(to_plot_rank),]
color_palette <- colorRampPalette(c("white", "red"))(100)
pheatmap::pheatmap(normalized_to_plot,cluster_rows = F,color = color_palette,cluster_cols = F)


