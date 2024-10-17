rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)

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
  Pancreas="panceas",
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
  negative_te <- negative[[tissue_public]]
  positive_te <- positive[[tissue_public]]
  annotation_te <- data.frame(te=c(positive_te,negative_te),condition=c(rep("Positive",length(positive_te)),rep("Negative",length(negative_te))))
  rownames(annotation_te) <- annotation_te$te
  annotation_te <- annotation_te[,-1,drop=F]
  to_plot <- data.frame()
  for(tissue in tissues){
    df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_TE.csv"),row.names = 1)
    positive_te_df <- df[which(rownames(df) %in% positive_te),]
    negative_te_df <- df[which(rownames(df) %in% negative_te),]
    search_table <- read.csv(paste0("data/samples/all/RNA_search_table.csv"))
    search_table <-search_table[which(search_table$tissue==tissue),]
    to_plot <- rbind(positive_te_df[,c(1:nrow(search_table)),drop=F],negative_te_df[,c(1:nrow(search_table)),drop=F])
    annotation_sample <- search_table[,c(3,5)]
    rownames(annotation_sample) <- annotation_sample$sample_name
    annotation_sample <- annotation_sample[,-1,drop=F]
    }
  max_abs <- max(abs(to_plot))
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

  heatmap_output <- pheatmap::pheatmap(to_plot,cluster_rows = F,cluster_cols = F,
                                       annotation_row = annotation_te,annotation_names_row = F,
                                       fontsize_col = 14,fontsize = 7,color = colors,breaks = breaks)
  print(heatmap_output)
}
# for(tissue_public in tissues_public){
#   public_TE_in_our_own_data(tissue_public)
# }
