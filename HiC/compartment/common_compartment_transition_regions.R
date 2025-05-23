rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
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
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus")
A_B_summary <- data.frame()
B_A_summary <- data.frame()
for(tissue in tissues){
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_50000.csv"))
  compartment$label <- paste0(compartment$chr,":",compartment$start,"-",compartment$end)
  compartment$tissue <- tissue_label_change(tissue)
  A_B <- compartment[which(compartment$condition=="A-B"),c("label","tissue")]
  B_A <- compartment[which(compartment$condition=="B-A"),c("label","tissue")]
  A_B_summary <- rbind(A_B_summary,A_B)  
  B_A_summary <- rbind(B_A_summary,B_A)
}
A_B_summary_count <- A_B_summary %>%   
  count(label)
A_B_summary_tissue <- A_B_summary %>%   
  group_by(label) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
A_B_summary_count <- merge(A_B_summary_count,A_B_summary_tissue,by="label")

B_A_summary_count <- B_A_summary %>%   
  count(label)
B_A_summary_tissue <- B_A_summary %>%   
  group_by(label) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
B_A_summary_count <- merge(B_A_summary_count,B_A_summary_tissue,by="label")

to_plot <- as.data.frame(table(A_B_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "Distribution of tissues number in common compartment A to B transition",  
    x = NULL,  
    y = "Count"  
  ) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme_minimal() +  
  xlab(NULL) +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 

to_plot <- as.data.frame(table(B_A_summary_count$n))
ggplot(data = to_plot, aes(x = Var1, y = Freq)) +  
  geom_bar(stat = "identity") +  
  labs(  
    title = "Distribution of tissues number in common compartment B to A transition",  
    x = NULL,  
    y = "Count"  
  ) +  
  geom_text(  
    aes(label = Freq),   
    vjust = 0,   
    size = 3.5  # Adjust text size as needed  
  ) +
  theme_minimal() +  
  xlab(NULL) +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  ) 
