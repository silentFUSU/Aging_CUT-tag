rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
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
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus","skin","muscle","cecum","ileum","pancreas","spleen")

summary <- data.frame()
resolution <- "25000_optimal_parameter"
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/diff_interaction_within_",resolution,"_loop.csv"))
  df <- as.data.frame(table(df$Significant))
  df$percent <- df$Freq/sum(df$Freq) *100
  df$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,df)
}
result <- summary %>%
  filter(Var1 %in% c("Down", "Up")) %>%
  group_by(tissue) %>%
  summarise(total_percent = sum(percent))
for(tissue in tissues){
  if(tissue_label_change(tissue) %in% result$tissue){
  }else{
    t_result <- data.frame(tissue=tissue_label_change(tissue),total_percent=0)
    result <- rbind(result,t_result)
  }
}
result <- result[order(result$total_percent),]
summary$tissue <- factor(summary$tissue,levels = result$tissue)
summary$Var1 <- factor(summary$Var1,levels = rev(c("Stable","Up","Down")))
color <- read.csv("data/samples/7_distinct_color.txt",header = F)
color <- setNames(c("gray","#e64b35","#3c5488"),c("Stable","Up","Down"))
p <- ggplot(summary, aes(x = percent, y = tissue, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Percent", y = "Tissue",fill="condition",title = "Loop") +
  scale_fill_manual(values = color)+
  theme_minimal()+
  theme(
    axis.text.x = element_text( size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal" 
  )+
  guides(fill = guide_legend(reverse = TRUE))             
ggsave("result/figures/loop_change_percentage.pdf",p,width = 6,height = 10)
