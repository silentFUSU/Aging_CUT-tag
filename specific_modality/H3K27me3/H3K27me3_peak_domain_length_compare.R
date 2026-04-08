rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
library(ggsignif)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
tissue_summary <- data.frame()
for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/edd_peaks_fdr05.bed"))
  domain$length <- domain$V3 - domain$V2 + 1
  domain_length <- mean(domain$length)
  
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_old_merge-W5000-G10000-E100.bed"))
  peak$length <- peak$V3 - peak$V2 + 1
  peak_length <- mean(peak$length)
  
  t_tissue_summary <- data.frame(tissue=tissue_label_change(tissue),peak_length=peak_length,domain_length=domain_length)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}
to_plot <- tissue_summary
to_plot <- reshape2::melt(to_plot)

breaks = seq(3,8,1)
labels <- ifelse(10**breaks >= 1e9, paste0(round(10**breaks/1e9, 1), "Gb"),
                 ifelse(10**breaks >= 1e6, paste0(round(10**breaks/1e6, 1), "Mb"), 
                        ifelse(10**breaks >= 1e3, paste0(round(10**breaks/1e3, 1), "Kb"),
                               as.character(10**breaks))))

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$tissue)))
p <- ggplot(to_plot, aes(x = variable , y = log10(value))) +
  geom_violin(fill = "lightgray", color = NA, alpha = 0.5) +
  geom_boxplot(outliers = F,width=0.1) +
  geom_jitter(aes(color = tissue), width = 0.2, size = 2) +
  scale_color_manual(values = color) +
  scale_y_continuous(breaks=breaks,
                     labels=labels)+
  labs(y = "Length") +
  theme_bw()+
  theme(
    axis.title.x = element_blank()
  )
ggsave("result/Sup_figures/H3K27me3_peak_domain_length_compare.pdf",p,width = 8,height = 8)
