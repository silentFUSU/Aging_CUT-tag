rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)  

options(scipen = 999)  
tissue <- "mammarygland"
resolution <- "50000"
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
check_row <- function(row) {  
  all(row[-1] == row[-1][1])  
}  


tissues <- c("lung","CB","liver","brain","kidney","colon","bonemarrow","thymus","heart","stomach","Hip","mammarygland")
to_plot <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(sample in search_table$sample_name){
    if(sample=="WJH-Liver-103"){
      df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1_recorrect.txt"))
    }else{
      df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    }
    # df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
    # df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X")))),]
    df <- df[,c(1,6)]
    colnames(df)[2] <- sample
    if(nrow(to_plot)==0){
      to_plot <- df
    }else{
      to_plot <- merge(to_plot,df,by="V1")
    }
  }
}
rownames(to_plot) <- to_plot$V1
to_plot <- to_plot[,-1]
pca <- prcomp(t(to_plot))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
search_table <- read.csv("data/samples/all/HiC_search_table.csv")
to_plot <- merge(to_plot,search_table,by="sample_name")
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
colours <- read.table("data/samples/30_distinct_color.txt")
colours <- setNames(colours$V1,tissues)
to_plot$age <- factor(to_plot$age,levels=c("3M","24M"))
ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue,shape=age)) + 
  scale_color_manual(values = colours) +
  geom_point(size=5) +
  theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  ggtitle("HiC")
