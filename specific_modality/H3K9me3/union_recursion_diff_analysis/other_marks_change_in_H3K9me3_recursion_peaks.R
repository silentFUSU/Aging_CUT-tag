rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
tissue <- "lung"
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
logFC_heatmap <- function(tissue,antibody){
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  decrease_H3K9me3 <- H3K9me3[which(H3K9me3$Significant=="Down"),c("Geneid","LogFC.old.young")]
  colnames(decrease_H3K9me3)[2] <- "H3K9me3"
  antibodys <- c("H3K27ac","H3K4me3","H3K4me1","H3K36me3","H3K27me3")
  for(antibody in antibodys){
    histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_H3K9me3_peaks_remove_batch_effect.csv"))
    histone <- histone[which(histone$Geneid %in% decrease_H3K9me3$Geneid), c("Geneid","LogFC.old.young")]
    colnames(histone)[2] <- antibody
    decrease_H3K9me3 <- merge(decrease_H3K9me3,histone,by="Geneid")
  }
  rna <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_in_H3K9me3_peaks.csv"))
  colnames(rna)[1] <- "Geneid"
  rna <- rna[which(rna$Geneid %in% decrease_H3K9me3$Geneid),c("Geneid","logFC")]
  colnames(rna)[2] <- "RNA"
  decrease_H3K9me3 <- merge(decrease_H3K9me3,rna,by="Geneid")
  rownames(decrease_H3K9me3) <- decrease_H3K9me3$Geneid
  breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
  decrease_H3K9me3 <- decrease_H3K9me3[order(decrease_H3K9me3$H3K9me3),]
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
  pheatmap::pheatmap(decrease_H3K9me3[,-1],breaks = breaks,cluster_cols = F,show_rownames = F,main = tissue_label_change(tissue),color = color_palette, )
}

antibody <- "H3K4me3"
logFC_heatmap_all_tissues <- function(antibody){
  H3K9me3 <- read.csv(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv"))
  H3K9me3 <- H3K9me3[order(H3K9me3$cluster),]
  tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
               "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","MEF")
  summary <- data.frame()
  for(tissue in tissues){
    if(antibody == "H3K9me3"){
      histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
    }else if(antibody=="ATAC"){
      histone <- read.csv(paste0("data/samples/",antibody,"/",tissue,"/",antibody,"/",antibody,"_diff_in_H3K9me3_peaks_remove_batch_effect.csv"))
    }else{
      histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_H3K9me3_peaks_remove_batch_effect.csv"))
    }
    histone <- histone[which(histone$Geneid %in% H3K9me3$X), c("Geneid","LogFC.old.young")]
    colnames(histone)[2] <- tissue_label_change(tissue)
    if(nrow(summary)==0){
      summary <- histone
    }else{
      summary <- merge(summary,histone,by="Geneid",all=T)
    }
  }
  summary$Geneid <- factor(summary$Geneid,levels = H3K9me3$X)
  summary <- summary[order(summary$Geneid),]
  rownames(summary) <- summary$Geneid
  breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
  annotation <- H3K9me3
  rownames(annotation) <- annotation$X
  annotation <- annotation[,-1,drop=F]
  annotation$cluster <- as.character(annotation$cluster)
  # tissues_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung","Mammary Gland",
  #                    "Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
  tissues_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Mef","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary","Colon","Stomach","Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
  # summary <- summary[,-1]
  summary <- summary[,tissues_order]
  annotation_color <- list(cluster=setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),c(1:4)))
  pheatmap::pheatmap(summary,breaks = breaks,cluster_cols = F,show_rownames = F,annotation_row=annotation,annotation_colors = annotation_color,color = color_palette,cluster_rows = F,main = antibody,legend = F,annotation_legend = F)
}
