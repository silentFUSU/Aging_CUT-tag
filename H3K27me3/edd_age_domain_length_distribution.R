rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggalluvial)  
library(data.table)

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
tissues <- sort(c("brain","CB","kidney","liver","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","lung"))
summary <- data.frame()
for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/edd_peaks_fdr05.bed"))
  t_summary <- data.frame(length=(domain$V3 - domain$V2 + 1), tissue= tissue_label_change(tissue))
  summary <- rbind(summary,t_summary)
}

color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,sort(unique(summary$tissue)))
ggplot(summary, aes(x = log10(length), fill = tissue)) +
  geom_histogram(binwidth = .1, alpha = 0.4, position = 'identity') +
  theme_minimal() +
  scale_fill_manual(values = color)+
  labs(title = "Distribution of Domain Length",
       x = "log10(Length)",
       y = "Count") +
  theme(legend.title = element_blank())+
  annotate("text", x = log10(200000), y = 40, label = "Length = 200Kb", color = "red", angle = 90, vjust = -0.3) +
  geom_vline(xintercept = log10(200000), linetype = "dashed", color = "red") 
result <- summary %>%
  group_by(tissue) %>%
  summarize(total_length = sum(length))

genome <- read.table("~/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes")
genome <- genome[which(genome$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
genome <- sum(genome$V2)
result$percentage <- result$total_length/genome *100
