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
resolution <- "50000"
summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
  df <- as.data.frame(table(df$condition))
  df$percent <- df$Freq/sum(df$Freq) *100
  df$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,df)
}
summary$Var1 <- factor(summary$Var1,levels = rev(c("A-A","B-B","A-B","B-A")))
result <- summary %>%
  filter(Var1 %in% c("A-B", "B-A")) %>%
  group_by(tissue) %>%
  summarise(total_percent = sum(percent))
result <- result[order(result$total_percent),]
color <- read.csv("data/samples/7_distinct_color.txt",header = F)
color <- setNames(color$V1,c("A-A","B-B","A-B","B-A"))
summary$tissue <- factor(summary$tissue,levels = result$tissue)
ggplot(summary, aes(x = percent, y = tissue, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Percent", y = "Tissue",fill="condition",title = "Compartment") +
  scale_fill_manual(values = color)+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),  
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom", 
    legend.direction = "horizontal" 
  )+
  guides(fill = guide_legend(reverse = TRUE))             

to_plot <- summary[which(summary$Var1 %in% c("A-B","B-A")),]
to_plot$percent[which(to_plot$Var1=="B-A")] <- -to_plot$percent[which(to_plot$Var1=="B-A")]
color <- setNames(c("#e64b35","#3c5488"),c("A-B","B-A"))
p <- ggplot(to_plot, aes(x = percent, y = tissue, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Percent", y = "Tissue",fill="condition",title = "Compartment") +
  scale_fill_manual(values = color)+
  theme_bw()+
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
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_continuous(limits = c(-2.5, 2.5),labels = function(x) format(abs(x)))
ggsave("result/figures/compartment_change_percentage.pdf",p,width = 6,height = 10)  

