rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(ggrepel)
library(Seurat)
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Frontal Cortex"
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
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return("10kb")
  }else{
    return("10kb")
  }
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
per_markers_all_tissue_pca <- function(antibody){
  if(antibody == "ATAC"){
    tab = read.delim(paste0("data/samples/ATAC/all/",antibody,"/merge-",bin_size(antibody),"_bins.counts"),row.names = 1,skip=1)    
  }else{
    tab = read.delim(paste0("data/samples/all/",antibody,"/merge-",bin_size(antibody),"_bins.counts"),row.names = 1,skip=1)    
  }
  tab <- tab[-which(tab$Chr=="chrY"),]
  counts <- tab[,c(6:ncol(tab))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
  if(antibody=="ATAC"){
    search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
  }else{
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    search_table <- search_table[which(search_table$antibody==antibody),]
    print(paste0(antibody," previous ",nrow(search_table)))
  }
  
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  search_table$tissue <- sapply(search_table$tissue,tissue_label_change)

  
  counts <- counts[,which(colnames(counts) %in% search_table$sample_name)]
  
  search_table$age[which(search_table$age=="3m")] <- "young"
  search_table$age[which(search_table$age=="24m")] <- "old"
  
  df_pca <- as.data.frame(edgeR::cpm(counts))
  variances <- apply(df_pca, 1, var)
  high_var_features <- names(sort(variances, decreasing = TRUE))[1:20000]
  selected_df_pca <- df_pca[high_var_features, ]
  pca_result <- prcomp(t(selected_df_pca), center = TRUE, scale. = TRUE)
  selected_pcs <- pca_result$x[, 1:30]
  umap_result <- umap::umap(selected_pcs,random_state=42)
  umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], names = colnames(df_pca))
  umap_df <- merge(umap_df,search_table[,c("sample_name","tissue","age")],by.x="names",by.y="sample_name")
  umap_df$tissue <- sapply(umap_df$tissue,tissue_label_change)
  umap_df$age <- factor(umap_df$age,levels=c("young","old"))
  

  colours <- read.table("data/samples/30_distinct_color.txt")
  colours <- setNames(colours$V1,sort(unique(search_table$tissue)))
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,color=tissue,shape=age)) + 
    scale_color_manual(values = colours) +
    geom_point(size=4) +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(antibody)
  
  return(p)
}
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","ATAC","H3K27ac","H3K4me1","H3K4me3")
p_list <- list()
for(i in c(1:length(antibodys))){
  p_list[[i]] <- per_markers_all_tissue_pca(antibodys[i])
}
combined_plot <- plot_a_list(p_list,no_of_rows = 2,no_of_cols = 4)
ggsave("result/Sup_figures/per_markers_all_tissues_UMAP_pca_select.pdf",combined_plot,width = 25,height = 10)
