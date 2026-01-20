rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyr)
library(dplyr)
library(tidyverse)
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

load("data/samples/GRN/grn_union_tissue.rdata")
load("data/samples/GRN/grn_union_skin.rdata")
grn_tissue[["skin"]] <- grn_union
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_list <- list()
for(tissue in tissues){
  grn <- grn_tissue[[tissue]]
  tissue_list[[tissue_label_change(tissue)]] <- unique(paste0(grn$peak,"-",grn$gene))
  }

n <- length(tissue_list)
similarity_matrix <- matrix(0, n, n)
rownames(similarity_matrix) <- names(tissue_list)
colnames(similarity_matrix) <- names(tissue_list)

for (i in 1:n) {
  for (j in 1:n) {
      intersection <- length(intersect(tissue_list[[i]], tissue_list[[j]]))
      union <- length(union(tissue_list[[i]], tissue_list[[j]]))
      similarity <- intersection / union
      similarity_matrix[i, j] <- similarity
  }
}

breaks <- c(seq(0.3, 1, length.out = 20))
color_palette <- colorRampPalette(c("white", "#8b0000"))(20)  
pheatmap::pheatmap(similarity_matrix,breaks = breaks,color = color_palette,filename = "result/Sup_figures/GRN_Cre_GENE_jaccard_similarity.pdf",width = 7,height = 6)


