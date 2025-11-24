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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
tissue_summary <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  df$length <- df$V3 - df$V2 +1
  df$tissue <- tissue_label_change(tissue)
  df <- df[,c("tissue","length")]
  tissue_summary <- rbind(tissue_summary,df)
  }
to_plot <- tissue_summary
to_plot$tissue <- factor(to_plot$tissue,levels=c("Cortex","Mammary Gland","Ovary","Muscle","Aorta","Hippocampus","BAT","Bladder","Lung","Thymus",
                                                 "Uterus","Cerebellum","Heart","Bone Marrow","Liver","iWAT","Tongue","Jejunum","Skin","Stomach",
                                                 "Testis","Spleen","Ileum","Colon","Kidney","Cecum","Pancreas"))

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))
ggplot(to_plot, aes(x = tissue, y = length,fill=tissue)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = color)+
  labs(title = "Boxplot of Length by Tissue",
       x = "Tissue",
       y = "Length") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))

