rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
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

df <- read.csv("data/samples/all/HiC_Quality_control.csv")
to_plot <- df[,c(1:4)]
to_plot$Age[which(to_plot$Age=="3M")] <-"Young"
to_plot$Age[which(to_plot$Age=="24M")] <- "Old"
to_plot <- to_plot %>%  
  mutate(TissueName = sapply(TissueName, tissue_label_change))  
to_plot$cis <- df$VP_cis/df$VP_unique * 100
to_plot$longrange_cis <- df$VP_cis...20K/df$VP_unique *100

to_plot$Age <- factor(to_plot$Age, levels=c("Young","Old"))
to_plot <- to_plot %>%
  arrange(TissueName,Age)
to_plot$SampleID <- paste0(to_plot$Age,"-",to_plot$TissueName,"-",to_plot$Mouse.ID)
to_plot$SampleID <- factor(to_plot$SampleID, levels=to_plot$SampleID)
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(to_plot$TissueName)))
ggplot(to_plot, aes(x = SampleID, y = cis, fill = TissueName, alpha = Age)) +  
  geom_col() +  
  scale_alpha_manual(values = c(Young = 1, Old = 0.6)) +
  labs(title = "Valid Cis / Valid pairs unique",  
       x = "Tissue",  
       y = "Cis%",  
       fill = "Tissue",   
       alpha = "Age Group") +  
  scale_fill_manual(values = color) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 10),legend.title = element_blank())

ggplot(to_plot, aes(x = SampleID, y = longrange_cis, fill = TissueName, alpha = Age)) +  
  geom_col() +  
  scale_alpha_manual(values = c(Young = 1, Old = 0.6)) +
  labs(title = "Valid Long range Cis / Valid pairs unique",  
       x = "Tissue",  
       y = "Long range Cis%",  
       fill = "Tissue",   
       alpha = "Age Group") +  
  scale_fill_manual(values = color) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 10),legend.title = element_blank())
