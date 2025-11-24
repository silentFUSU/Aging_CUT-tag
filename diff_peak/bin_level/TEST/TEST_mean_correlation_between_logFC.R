rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(corrplot)  
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

antibodys <- c("H3K9me3","H3K27me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
matrix_size <- length(tissues)
correlation_matrix <- matrix(NA,nrow=matrix_size,ncol=matrix_size)
dimnames(correlation_matrix) <- list(tissues, tissues)  
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
      tissue_label <- "Mammarygland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

correlation_matrix <- as.data.frame(correlation_matrix)
bin_size <- function(antibody){
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    return("10kb")
  }else{
    return("1kb")
  }
}

for(i in c(1:(length(tissues)-1))){
  tissue_row <- tissues[i]
  for(j in c((i+1):length(tissues))){
    tissue_col <- tissues[j]
    print(paste0(tissue_row,"-",tissue_col))
    cor <- 0
    for(antibody in antibodys){
      tissue_row_df <- read.csv(paste0("data/samples/",tissue_row,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff.csv"))
      tissue_col_df <- read.csv(paste0("data/samples/",tissue_col,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff.csv"))
      df <- merge(tissue_row_df[,c("Geneid","LogFC.old.young")],tissue_col_df[,c("Geneid","LogFC.old.young")],by="Geneid")
      t_cor <- cor(df[,c(2,3)])
      t_cor <- t_cor[1,2]
      cor <- cor + t_cor
    }
    cor <- cor / length(antibodys)
    correlation_matrix[i,j] <- cor
  }
}
# saveRDS(correlation_matrix,"data/samples/all/all_logFC_correlation.rds")
correlation_matrix <- readRDS("data/samples/all/all_logFC_correlation.rds")
for(i in c(1:nrow(correlation_matrix))){
  correlation_matrix[i,i] <- 1
}
for(i in c(2:nrow(correlation_matrix))){
  for(j in c(1:(i-1))){
    correlation_matrix[i,j] <- correlation_matrix[j,i]
  }
}
rownames(correlation_matrix) <- sapply(rownames(correlation_matrix), tissue_label_change)  
colnames(correlation_matrix) <- sapply(colnames(correlation_matrix), tissue_label_change)  

pheatmap::pheatmap(correlation_matrix,breaks = seq(-0.4, 0.4, length.out = 101))
