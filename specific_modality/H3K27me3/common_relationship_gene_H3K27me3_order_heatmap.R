rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
summary <- data.frame()
annotation <- read.csv("data/samples/RNA/H3K27me3_decreased_RNA_increased_H3K27me3_union_kmeans.csv")
for(tissue in tissues){
  RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  colnames(RNA)[1] <- "Geneid"
  RNA <- RNA[which(RNA$Geneid %in% annotation$X),c("Geneid","logFC")]
  colnames(RNA)[2] <- tissue_label_change(tissue)
  if(nrow(summary)==0){
    summary <- RNA
  }else{
    summary <- merge(summary,RNA,by="Geneid",all=T)
  }
}
colnames(annotation)[1] <- "Geneid"
summary <- merge(summary,annotation,by="Geneid")
rownames(summary) <- summary$Geneid
summary <- summary[order(summary$cluster),]
summary <- summary[,-c(1,ncol(summary))]
rownames(annotation) <- annotation$Geneid
annotation <- annotation[,-1,drop=F]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))  
annotation$cluster <- as.character(annotation$cluster)
# tissues_order <- c("Ovary","Thymus","Uterus","Muscle","Bladder","Liver","Kidney","Mammary Gland","Ileum","Cecum","Cortex","Hippocampus","Bone Marrow","IWAT","Spleen","Jejunum","Stomach","Pancreas","Testis","Heart","Tongue","Aorta","BAT","Lung","Colon","Cerebellum","Skin")
tissues_order <- c("Ovary","Uterus","Muscle","Thymus","Kidney","Liver","Mammary Gland","Bone Marrow","IWAT","Jejunum","Spleen","Cecum","Ileum","Cerebellum","Skin","Cortex","Hippocampus","Bladder","Testis",
                   "Pancreas","Stomach","Aorta","BAT","Colon","Lung","Heart","Tongue")
summary <- summary[,tissues_order]
pheatmap::pheatmap(summary,cluster_rows = F,cluster_cols = F,show_rownames = F,na_col = "Grey",breaks = breaks, annotation_row = annotation, color = color_palette)

