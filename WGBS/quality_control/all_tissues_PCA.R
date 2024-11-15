rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(Polychrome)
tissues <- c("mammarygland","lung","liver","kidney","ileum","Hip","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis","heart","stomach","muscle","bladder","aorta","tongue","spleen","pancreas","brain","cecum","uterus","iWAT")
bin_size <- "1kb"

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
    }
  }
  return(tissue_label)
}

all_tissues_PCA <- function(tissues,bin_size){
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    t_df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"))
    t_df <- reshape2::dcast(t_df,label ~ sample, value.var = "percent")
    if(i == 1){
      df <- t_df
    }else{
      df <- merge(df,t_df,by="label")
    }
  }
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  rownames(df) <- df$label
  df <- df[,-1]
  df <- na.omit(df) 
  pca <- prcomp(t(df))
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  to_plot <- merge(to_plot,search_table,by="sample_name")
  to_plot$age[which(to_plot$age == "3M")] <- "young"
  to_plot$age[which(to_plot$age == "24M")] <- "old"
  to_plot$age <- factor(to_plot$age,levels = c("young","old"))
  for(j in c(1:nrow(to_plot))){
    to_plot$tissue[j] <- tissue_label_change(to_plot$tissue[j])
  }
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  color_tissues <- sapply(sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                                 "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")), tissue_label_change) 
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- setNames(color$V1,color_tissues)
  
  p<- ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
    geom_point(size=5) +theme_bw()+
    scale_color_manual(values = color) +
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = mouse_ID, color = tissue),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    )+ggtitle(paste0("compress to ",bin_size," bin"))
  return(p)
}
bin_size <- "1kb"
p <- all_tissues_PCA(tissues,bin_size)
p <- readRDS("p.rds")
