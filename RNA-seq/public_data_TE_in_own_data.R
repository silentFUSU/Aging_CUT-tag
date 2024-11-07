rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(dplyr)
library(stringr)
conditions <- c("Negative","Positive")
tissues_public <- c("BAT","Bone","Brain","GAT","Heart","Kidney","Limb_Muscle","Liver","Lung",
             "MAT","Marrow","Pancreas","SCAT","Skin","Small_Intestine","Spleen","WBC")
positive <- readRDS("data/public_data/GSE132040/TE/positive_TE_list_pvalue005.rds")
negative <- readRDS("data/public_data/GSE132040/TE/negative_TE_list_pvalue005.rds")

dict <- list(
  BAT="BAT",
  Heart="heart",
  Kidney="kidney",
  Limb_Muscle="muscle",
  Liver="liver",
  Lung="lung",
  Marrow="bonemarrow",
  Pancreas="pancreas",
  Skin="skin",
  Brain=c("FC","CB","Hip"),
  Spleen="spleen",
  Small_Intestine=c("ileum","jejunum")
)
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if (tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

public_TE_in_our_own_data <- function(tissue_public){
  tissues <- dict[[tissue_public]]
  if(length(tissues)==0){
    return()
  }
  top <- 10
  negative_te <- negative[[tissue_public]]
  negative_te <- negative_te[order(negative_te$p_value),]
  positive_te <- positive[[tissue_public]]
  positive_te <- positive_te[order(positive_te$p_value),]
  
  annotation_te <- data.frame(te=c(positive_te$TE[1:min(top,nrow(positive_te))],negative_te$TE[1:min(top,nrow(negative_te))]),condition=c(rep("Positive",min(top,nrow(positive_te))),rep("Negative",min(top,nrow(negative_te)))))
  rownames(annotation_te) <- annotation_te$te
  annotation_te <- annotation_te[,-1,drop=F]
  to_plot <- data.frame()
  for(tissue in tissues){
    df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_TE.csv"),row.names = 1)
    colnames(df)[which(colnames(df)=="logFC")] <- tissue_label_change(tissue)
    positive_te_df <- df[which(rownames(df) %in% rownames(annotation_te)[which(annotation_te$condition=="Positive")]),]
    negative_te_df <- df[which(rownames(df) %in% rownames(annotation_te)[which(annotation_te$condition=="Negative")]),]
    # search_table <- read.csv(paste0("data/samples/all/RNA_search_table.csv"))
    # search_table <-search_table[which(search_table$tissue==tissue_label_change(tissue)),]
    t_to_plot <- rbind(positive_te_df[,tissue_label_change(tissue),drop=F],negative_te_df[,tissue_label_change(tissue),drop=F])
    if(nrow(to_plot)==0){
      to_plot <- t_to_plot
    }else{
      to_plot <- merge(to_plot,t_to_plot,by="row.names")
      rownames(to_plot) <- to_plot$Row.names
      to_plot <- to_plot[,-1]
    }
    # annotation_sample <- search_table[,c(3,5)]
    # rownames(annotation_sample) <- annotation_sample$sample_name
    # annotation_sample <- annotation_sample[,-1,drop=F]
  }
  to_plot$rownames <- rownames(to_plot)
  to_plot$rownames <- factor(to_plot$rownames,levels=rownames(annotation_te))
  to_plot <-  to_plot[order(to_plot$rownames),]
  to_plot <- to_plot[,-ncol(to_plot),drop=F]
  max_abs <- max(abs(to_plot))
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  # colors <- colorRampPalette(c("#436CA7", "#FFFCBF", "white","#FFFCBF","#C5392B"))(length(breaks) - 1)
  heatmap_output <- pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,
                                       annotation_row = annotation_te,annotation_names_row = F,
                                       fontsize_col = 14,fontsize_row = 10,fontsize = 7,breaks = breaks)
  
  # heatmap_output <- pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,
  #                                      annotation_row = annotation_te,annotation_names_row = F,
  #                                      fontsize_col = 14,fontsize = 7)
  print(heatmap_output)
}
for(tissue_public in tissues_public){
  public_TE_in_our_own_data(tissue_public)
}
