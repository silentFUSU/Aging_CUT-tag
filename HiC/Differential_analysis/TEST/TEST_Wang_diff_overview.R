rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
resolution <- "200000"
tissues <- c("brain","CB","liver","lung","kidney","colon","heart","bonemarrow","stomach","thymus","mammarygland","Hip")
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
to_plot <- data.frame()
for(tissue in tissues){
  re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
  re$Significant <- "Stable"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
  re <- re[which(abs(re$V2-re$V3)>5),]
  t_to_plot <- as.data.frame(table(re$Significant))
  t_to_plot <- t_to_plot[which(t_to_plot$Var1 %in% c("Down","Up")),]
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq)*100
  t_to_plot$tissue <- tissue_label_change(tissue)
  to_plot <- rbind(to_plot,t_to_plot)
}

ggplot(to_plot, aes(x = tissue, y = Freq, fill = Var1)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ggtitle("Differential interactions")

ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Percent")+
  ggtitle("Differential interactions")

## chromosome propotion
to_plot <- data.frame()
for(tissue in tissues){
  re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
  re$Significant <- "Stable"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
  re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
  re <- re[which(abs(re$V2-re$V3)>5),]
  re$V1 <- paste0("chr",re$V1)
  re$V1[which(re$V1=="chr20")] <- "chrX"
  re$V1[which(re$V1=="chr21")] <- "chrY"
  re <- re[which(re$Significant != "Stable"),]
  increase <- re[which(re$Significant=="Up"),]
  increase <- as.data.frame(table(increase$V1))
  increase$increase_percent <- increase$Freq/sum(increase$Freq)*100  
  
  decrease <- re[which(re$Significant=="Down"),]
  decrease <- as.data.frame(table(decrease$V1))
  decrease$decrease_percent <- decrease$Freq/sum(decrease$Freq)*100
  
  t_to_plot <- merge(increase[,c("Var1","increase_percent")], decrease[,c("Var1","decrease_percent")],by="Var1")
  t_to_plot$tissue <- tissue_label_change(tissue)
  
  to_plot <- rbind(to_plot,t_to_plot)
}
to_plot$Var1 <- factor(to_plot$Var1, levels=paste0("chr",c(1:19,"X","Y")))
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,paste0("chr",c(1:19,"X","Y")))
ggplot(to_plot, aes(x = tissue, y = increase_percent, fill = Var1)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_manual(values=color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ggtitle("Increased interactions")

ggplot(to_plot, aes(x = tissue, y = decrease_percent, fill = Var1)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_manual(values=color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ggtitle("Decreased interactions")
