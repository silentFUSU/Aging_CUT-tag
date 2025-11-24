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
annotation <- read.csv("data/samples/all/H3K27me3/recursion_bin_diff_table/kmeans_annotation.csv")
ref <- read.table("~/ref_data/mm10_10kb_bins.bed")
annotation <- merge(annotation,ref,by.x="X",by.y="V4")

result <- annotation %>%
  group_by(cluster, V1) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count) *100) %>%
  ungroup()

result$V1 <- factor(result$V1,levels = paste0("chr",c(1:19,"X","Y")))
color <- read.csv("data/samples/30_distinct_color.txt",header = F)
color <- setNames(color$V1,paste0("chr",c(1:19,"X","Y")))
result$cluster <- paste0("kmeans",result$cluster)
ggplot(result, aes(x = cluster, y = proportion, fill = V1)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color,name = "Chromosome") +
  labs(title = "Proportion of Chromosomes within each Cluster",
       x = "Cluster",
       y = "Proportion") +
  theme_minimal()

