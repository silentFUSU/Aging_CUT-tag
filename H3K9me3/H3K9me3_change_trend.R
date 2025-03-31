rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)
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
tissues <- c("aorta","bonemarrow","brain","mammarygland","heart","Hip","lung","iWAT","ovary","skin","uterus")
tissues <- c("testis","tongue")
df <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_bedtools_filtered.bed")
df <- df[-which(df$V1=="chrY"),]
antibody <- "H3K9me3"
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  # df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
  df <- df[which(df$Significant!="Stable"),]
  df <- df[,c("Significant","LogFC.old.young")]
  df$LogFC.old.young <- abs(df$LogFC.old.young)
  df$Significant <- factor(df$Significant,levels=c("Up","Down"))
  pvalue <- t.test(abs(df$LogFC.old.young[which(df$Significant=="Up")]), abs(df$LogFC.old.young[which(df$Significant=="Down")]))
  pvalue <- pvalue$p.value
  p <- ggplot(df, aes(x = Significant, y = LogFC.old.young, fill= Significant)) +  
    geom_violin(adjust = 2.5) +          
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  
    theme_minimal()+
    ggtitle(tissue_label_change(tissue))+
    theme(text = element_text(size = 20),legend.position = "none")+
    labs(x = NULL,y = "abs(log2(Fold Change))")+
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(pvalue, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red")
  
  print(p)
}
antibody <- "H3K9me3"
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
  search_table$age[which(search_table$age=="3m")] <- "young"
  search_table$age[which(search_table$age=="24m")] <- "old"
  samples <- paste(search_table$sample_name,search_table$age,search_table$mouse_ID,search_table$batch,sep = ".")
  df <- df[,c("Geneid",samples)]
  rownames(df) <- df$Geneid
  df <- df[,-1]
  df <- reshape2::melt(df)
  search_table$variable <- paste(search_table$sample_name,search_table$age,search_table$mouse_ID,search_table$batch,sep = ".")
  df <- merge(df,search_table,by="variable")
  p <- ggplot(df, aes(x = age, y = log2(value), fill= age)) +  
    geom_violin(adjust = 2.5) +          
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  
    theme_minimal()+
    ggtitle(tissue_label_change(tissue))+
    theme(text = element_text(size = 20),legend.position = "none")+
    labs(x = NULL,y = "abs(log2(Fold Change))")
  t.test(df$value[which(df$age=="young")], df$value[which(df$age=="old")])
  print(p)
}
