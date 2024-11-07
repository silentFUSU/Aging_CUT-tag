rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)

df <- read.csv("data/samples/all/WGBS_depth.csv")
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
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

df <- df %>%  
  mutate(TissueName = sapply(TissueName, tissue_label_change)) 
df$Depth..Million. <- df$Depth..Million.*2*150/1024/1024/1024
df$Age <- factor(df$Age, levels = c("3M","24M"))
ggplot(df,aes(x=TissueName,y=Depth..Million.,color = Age))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
  geom_text(aes(label = Mouse.ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
  scale_fill_brewer(palette="Set3")+
  ggtitle(paste0("WGBS sequencing depth"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Depth (G)")

ggplot(df,aes(x=TissueName,y=Depth..Million.,fill=TissueName))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.7)+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  theme_bw()+
  theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Depth (G)")+
  xlab("")+labs(fill = "", color = "")+  
  ggtitle(paste0("WGBS sequencing depth"))

mean(df$Depth..Million.)
