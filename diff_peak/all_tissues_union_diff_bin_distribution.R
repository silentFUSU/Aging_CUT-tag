rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
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
antibody <- "H3K27me3"
common_increase <- read.csv(paste0("data/samples/all/",antibody,"/common_increase_10kb_bins_after_remove_batch_effect.csv"))
common_decrease <- read.csv(paste0("data/samples/all/",antibody,"/common_decrease_10kb_bins_after_remove_batch_effect.csv"))
common_increase <- as.data.frame(table(common_increase$n))
common_decrease <- as.data.frame(table(common_decrease$n))
# common_decrease$Freq <- -common_decrease$Freq
common_increase$condition <- "common increase"
common_decrease$condition <- "common decrease"
common_increase$percent <- common_increase$Freq/sum(common_increase$Freq)*100
common_decrease$percent <- common_decrease$Freq/sum(common_decrease$Freq)*100
common_df <- rbind(common_increase,common_decrease)
common_df$condition <- factor(common_df$condition,levels=c("common increase","common decrease"))
ggplot() +  
  geom_line(data = common_df, aes(x = Var1, y = percent, color = condition, group = condition), size = 1) +  
  labs(  
    title = "Distribution of tissues number in common changed bins",  
    x = NULL,  
    y = "Percnet"  
  ) +  
  theme_minimal() +  
  xlab(NULL) +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) +  
  ylim(0,30)  


ggplot(common_df, aes(x = Var1, y = Freq, fill = condition)) +  
  geom_bar(stat = "identity") +  
  labs(title = paste0("Distribution of tissues number in common changed bins"), x = NULL, y = "Count") +  
  theme_minimal() +  
  xlab(NULL)+
  scale_y_continuous(labels = abs) +  
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold")
  ) +  
  ylim(-48000,48000)+
  scale_fill_manual(values = c("common decrease" = "skyblue",  "common increase"= "salmon"), name = NULL) 

common_decrease$Freq <- -common_decrease$Freq
ggplot(common_increase, aes(x = Var1, y = Freq, fill = condition)) +  
  geom_bar(stat = "identity",alpha=0.5) +  
  geom_bar(data = common_decrease, aes(x = Var1, y = Freq), stat = "identity", fill = "skyblue", alpha = 0.5) +  
  labs(  
    title = "Distribution of tissues number in common changed bins",  
    x = NULL,  
    y = "Count"  
  ) +  
  theme_minimal() +  
  xlab(NULL) +  
  scale_y_continuous(labels = abs) +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) +  
  ylim(0, 48000)  


