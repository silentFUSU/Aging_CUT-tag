rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
state <- "14"

model_plot <- function(tissues,state){
  model <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state,"_all_tissues/emissions_",state,".txt"))
  model <- model[,-1]
  model <- model[,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")]
  rownames(model) <- paste0("state",c(1:nrow(model)))
  color_palette <- colorRampPalette(c("white", "blue"))(50) 
  pheatmap::pheatmap(model,cluster_cols = F,cluster_rows = F,color = color_palette)
  
  
  percent_all <- data.frame(Var1 = as.character(),
                            Freq = as.numeric(),
                            tissue = as.character(),
                            age = as.character())
  ages <- c("young1","young2","old1","old2")
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    for(age in ages){
      df <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state,"_all_tissues/split_1k/",tissue,"_",age,"_",state,"_segments_1k.bed"),header = F)
      percent <- as.data.frame(table(df$V4))
      percent$Var1 <- factor(percent$Var1, levels = paste0("E",c(1:nrow(percent))))
      percent <- percent[order(percent$Var1),]
      percent$tissue <- tissue
      percent$age <- age
      percent_all <- rbind(percent_all,percent)
    }
  }
  percent_all_to_plot <-  percent_all %>%
    group_by(Var1) %>%
    summarise(Freq = sum(Freq))
  
  result <- percent_all_to_plot %>%  
    mutate(Percent = Freq / sum(Freq) * 100,  
           # 计算图表标签  
           Label = paste0(Var1, ": ", round(Percent, 1), "%")) 
  color <- read.table("data/samples/20_distinct_color.txt")
  color <- setNames(color$V1,result$Var1)
  ggplot(result, aes(x = "", y = Freq, fill = Var1)) +  
    geom_bar(stat = "identity", width = 1) +  
    coord_polar(theta = "y") +  
    theme_void() +  # 移除背景和坐标轴  
    scale_fill_manual(values = color, labels = result$Label ) +  
    theme(text = element_text(size = 20))+
    guides(fill = guide_legend(title = "State"))  
}


tissue <- "CB"
state <- "16"
young1 <- read.delim(paste0("result/all/ChromHMM/until_ovary/",state,"_until_ovary/split_1k/",tissue,"_young1_",state,"_segments_1k.bed"),header = F)
young2 <- read.delim(paste0("result/all/ChromHMM/until_ovary/",state,"_until_ovary/split_1k/",tissue,"_young2_",state,"_segments_1k.bed"),header = F)
old1 <- read.delim(paste0("result/all/ChromHMM/until_ovary/",state,"_until_ovary/split_1k/",tissue,"_old1_",state,"_segments_1k.bed"),header = F)
old2 <- read.delim(paste0("result/all/ChromHMM/until_ovary/",state,"_until_ovary/split_1k/",tissue,"_old2_",state,"_segments_1k.bed"),header = F)
young1 <- young1[which(young1$V4=="E16"),]
young2 <- young2[which(young2$V4=="E16"),]
old1 <- old1[which(old1$V4=="E16"),]
old2 <- old2[which(old2$V4=="E16"),]
young1$label <- paste0(young1$V1,"-",young1$V2,"-",young1$V3)
young2$label <- paste0(young2$V1,"-",young2$V2,"-",young2$V3)
old1$label <- paste0(old1$V1,"-",old1$V2,"-",old1$V3)
old2$label <- paste0(old2$V1,"-",old2$V2,"-",old2$V3)
young <- young1[which(young1$label %in% young2$label),]
old <- old1[which(old1$label %in% old2$label),]
merge <- young[which(young$label %in% old$label),]
