rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(maditr)
state_num <- 14
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
chromHMM_state_transfer_H3K27me3_H3K9me3_change <- function(tissue,state_num){
  H3K27me3_H3K9me3 <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant.bed"))
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  
  young$label <- paste0(young$V1,":",young$V2,"-",young$V3)
  old$label <- paste0(old$V1,":",old$V2,"-",old$V3)
  young <- young[which(young$label %in% old$label),]  
  old <- old[which(old$label %in% young$label),]
  colnames(young)[4] <- "young_state"
  colnames(old)[4] <- "old_state"
  state_change <- merge(young[,c(1:5)],old[,c(4:5)],by="label")
  state_change <- state_change[which(state_change$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
  state_change <- state_change[,c(-1)]
  state_change$condition <- paste0(state_change$young_state,"-",state_change$old_state)
  
  H3K27me3_H3K9me3$V2 <- H3K27me3_H3K9me3$V2 + 1
  H3K27me3_H3K9me3 <- as.data.table(H3K27me3_H3K9me3)
  setDT(H3K27me3_H3K9me3)  
  setkey(H3K27me3_H3K9me3,V1,V2,V3)
  
  state_change$V2 <- state_change$V2 + 1 
  state_change <- as.data.table(state_change)
  setDT(state_change) 
  setkey(state_change,V1,V2,V3)
  
  overlaps <- foverlaps(H3K27me3_H3K9me3,state_change, type = "any", nomatch = 0L)
  to_plot <- as.data.frame(table(overlaps$condition))
  to_plot <- to_plot[which(to_plot$Freq >100),]
  to_plot$percent <- to_plot$Freq/sum(to_plot$Freq)*100
  ggplot(to_plot, aes(x=1,y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)))
  }