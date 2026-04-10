rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)

search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
domain_list <- data.frame()
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K27me3"

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

summary <- data.frame()
for(tissue in tissues){
  domains <- read.table(paste0("data/samples/",tissue,"/",antibody,"/peaks/edd/edd_peaks_fdr05.bed"))
  domains$length <- domains$V3 - domains$V2 +1
  domains <- domains[which(domains$length >= 200000),]
  t_summary <- data.frame(tissue=tissue_label_change(tissue),count=nrow(domains),coverage=sum(domains$length))
  summary <- rbind(summary,t_summary)
}

to_plot <- summary
to_plot$label <- "domain"
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(to_plot$tissue))
p <- ggplot(to_plot, aes(x = label,y = count),) +
  geom_violin(fill = "lightgray", color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.2) +
  geom_jitter(aes(color = tissue), width = 0.2, size = 2) + 
  scale_color_manual(values = color) +
  labs(y = "count",x=NULL) +
  theme_bw()+
  theme(
    axis.title.x = element_blank(),  
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank()
  )+ylim(0,600)
p
ggsave("result/Sup_figures/H3K27me3_domain_number.pdf",p,width = 6,height = 8)
p<- ggplot(to_plot, aes(x = tissue, y = count,fill=tissue)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = color) +
  labs(
    x = NULL,
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

breaks = seq(7,10,1)
labels <- ifelse(10**breaks >= 1e9, paste0(round(10**breaks/1e9, 1), "Gb"),
                 ifelse(10**breaks >= 1e6, paste0(round(10**breaks/1e6, 1), "Mb"), 
                        ifelse(10**breaks >= 1e3, paste0(round(10**breaks/1e3, 1), "Kb"),
                               as.character(10**breaks))))

p <- ggplot(to_plot, aes(x = label , y = log10(coverage))) +
  geom_violin(fill = "lightgray", color = NA, alpha = 0.5) +
  geom_boxplot(outliers = F) +
  geom_jitter(aes(color = tissue), width = 0.2, size = 2) +
  scale_color_manual(values = color) +
  scale_y_continuous(breaks=breaks,
                     labels=labels)+
  labs(y = "coverage") +
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),  
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank()
  )
ggsave("result/Sup_figures/H3K27me3_domain_coverage.pdf",p,width = 6,height = 8)
ggplot(to_plot, aes(x = tissue, y = coverage,fill=tissue)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = color) +
  labs(
    x = NULL,
    y = "Coverage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

