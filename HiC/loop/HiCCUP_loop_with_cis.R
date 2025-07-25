rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
library(tidyr)
library(ggrepel)
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

resolution <- "25000_optimal_parameter"
quality_control <- read.csv("data/samples/all/HiC_Quality_control.csv")
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(sample in search_table$sample_name){
    df <- read.csv(paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/",sample,"_",resolution,"/",sample,"_",resolution,"_loop.csv"))    
    t_summary <- data.frame(SampleID=sample,tissue=tissue_label_change(tissue),counts=nrow(df))
    summary <- rbind(summary,t_summary)
    }
}

summary <- merge(summary,quality_control[,c("SampleID","Age","VP_cis.VP_unique")],by="SampleID")
summary$VP_cis.VP_unique <- as.numeric(gsub("%", "", summary$VP_cis.VP_unique))
summary$VP_cis.VP_unique <- as.numeric(summary$VP_cis.VP_unique)

ggplot(summary, aes(x =VP_cis.VP_unique , y =  counts , color = tissue,shape=Age)) +
  geom_point(size = 3) +  
  labs(
    title = "VP_cis/VP_unique vs. Counts",
    x = "VP_cis/VP_unique",
    y = "Loop Counts"
  ) + 
  theme_minimal()
