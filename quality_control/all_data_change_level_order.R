rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
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
  summary_per_antibody <- summary_per_antibody[order(summary_per_antibody$count,decreasing = TRUE),]
  summary_per_antibody$rank <- c(1:length(tissues))
  colnames(summary_per_antibody)[which(colnames(summary_per_antibody)=="rank")] <- antibody
  summary_per_antibody <- summary_per_antibody[,c(1,3)]
  if(nrow(summary)==0){
    summary <- summary_per_antibody
  }else{
    summary <- merge(summary,summary_per_antibody,by="tissue")
  }
}
row_sums <- rowSums(summary[, -1])  
sorted_indices <- order(row_sums)  
summary <- summary[sorted_indices, ]
rownames(summary) <- summary$tissue
breaks <- seq(1, 27, length.out = 27) 
color_palette <- colorRampPalette(c("red", "#ffe2e2","white"))(27)  
pheatmap::pheatmap(summary[,-1],cluster_rows = F,cluster_cols = F,breaks = breaks,color = color_palette)

##### gene_expression
summary_RNA <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  df <- df[which(df$Significant !="Stable"),]
  t_summary_RNA <- data.frame(tissue=tissue_label_change(tissue),count=nrow(df))  
  summary_RNA <- rbind(summary_RNA,t_summary_RNA)
}
summary_RNA <- summary_RNA[order(summary_RNA$count,decreasing = TRUE),]
summary_RNA$rank <- c(1:length(tissues))
colnames(summary_RNA)[3] <- "RNA"
summary_RNA <- summary_RNA[,c(1,3)]
summary <- merge(summary,summary_RNA,by="tissue")

row_sums <- rowSums(summary[, -1])  
sorted_indices <- order(row_sums)  
summary <- summary[sorted_indices, ]
rownames(summary) <- summary$tissue
breaks <- seq(1, 27, length.out = 27) 
color_palette <- colorRampPalette(c("red", "#ffe2e2","white"))(27)  
pheatmap::pheatmap(summary[,-1],cluster_rows = F,cluster_cols = F,breaks = breaks,color = color_palette)

##### DNA methylation
summary_DNA <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta0.txt"))
  t_summary_DNA <- data.frame(tissue=tissue_label_change(tissue),count=nrow(df))
  summary_DNA <- rbind(summary_DNA,t_summary_DNA)
}

summary_DNA <- summary_DNA[order(summary_DNA$count,decreasing = TRUE),]
summary_DNA$rank <- c(1:length(tissues))
colnames(summary_DNA)[3] <- "WGBS"
summary_DNA <- summary_DNA[,c(1,3)]
summary <- merge(summary,summary_DNA,by="tissue")

row_sums <- rowSums(summary[, -1])  
sorted_indices <- order(row_sums)  
summary <- summary[sorted_indices, ]
rownames(summary) <- summary$tissue
breaks <- seq(1, 27, length.out = 27) 
color_palette <- colorRampPalette(c("red", "#ffe2e2","white"))(27)  
pheatmap::pheatmap(summary[,-1],cluster_rows = F,cluster_cols = F,breaks = breaks,color = color_palette)
