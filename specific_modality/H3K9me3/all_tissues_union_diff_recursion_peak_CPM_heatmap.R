rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
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
common_increase <- read.csv(paste0("data/samples/all/",antibody,"/common_increase-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_increase <- common_increase[which((common_increase$end - common_increase$start +1) > 200000),]
common_decrease <- read.csv(paste0("data/samples/all/",antibody,"/common_decrease-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_decrease <- common_decrease[which((common_decrease$end - common_decrease$start +1) > 200000),]
regions <- c(common_increase$Geneid[which(common_increase$n>0)],common_decrease$Geneid[which(common_decrease$n>0)])
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
CPM_summary <- data.frame()

for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody=="H3K9me3" & search_table$age=="3m"),]
  colnames <- paste0(search_table$sample_name,".young.",search_table$mouse_ID,".",search_table$batch)
  df <- df[which(df$Geneid %in% regions),c("Geneid",colnames)]
  df$CPM <- rowSums(df[,-1])/length(colnames)
  df <- df[,c("Geneid","CPM")]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(CPM_summary) == 0){
    CPM_summary <- df
  }else{
    CPM_summary <- merge(CPM_summary,df,by="Geneid")    
  }
}
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
annotation <- annotation[order(annotation$cluster),]
rownames(CPM_summary) <- CPM_summary$Geneid
CPM_summary <- CPM_summary[,-1]
tissues_order <- c("Kidney","Muscle","Skin","Liver","Aorta","Testis","Cortex","Tongue","Uterus","Bladder","Stomach","Heart","Hippocampus","Cerebellum","BAT","Lung","Mammary Gland",
                   "Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
CPM_summary <- CPM_summary[annotation$X,tissues_order]
rownames(annotation) <- annotation$X
annotation <- annotation[,-1,drop=F]
CPM_summary <- log2(CPM_summary)
annotation$cluster <- as.character(annotation$cluster)
breaks <- c(seq(0, 5.9, length.out = 40), seq(6, 8.9, length.out = 20), seq(9, 15, length.out = 40))  
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
pheatmap::pheatmap(CPM_summary,cluster_rows = F,cluster_cols = F,annotation_row = annotation,breaks = breaks,show_rownames = F,color = color_palette)
pheatmap::pheatmap(CPM_summary,cluster_rows = F,cluster_cols = F,annotation_row = annotation,scale="row",show_rownames = F,color = color_palette)


