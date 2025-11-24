rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
options(scipen = 999)
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

calculate_overlap <- function(start1, end1, start2, end2) {
  max_overlap_start <- max(start1, start2)
  min_overlap_end <- min(end1, end2)
  overlap_length <- max(0, min_overlap_end - max_overlap_start)
  return(overlap_length)
}

PMD_HMD_region <- read.table("data/public_data/PMD_coordinates_mm10.bed")
PMD_HMD_region$V2 <- PMD_HMD_region$V2 + 1
PMD_HMD_region$label <- paste0("bin",c(1:nrow(PMD_HMD_region)))
PMD_HMD_region <- PMD_HMD_region[,c("V1","V2","V3","V5","label")]
PMD_HMD_region$V5[is.na(PMD_HMD_region$V5)] <- "other"
PMD_HMD_region <- as.data.table(PMD_HMD_region)
setDT(PMD_HMD_region)
setkey(PMD_HMD_region,V1,V2,V3)

H3K27me3_peaks <- read.csv("data/samples/all/H3K27me3/peaks_merged/kmeans_annotation_larger_20_tissues_RPKM.csv")
split_chr <- strsplit(as.character(H3K27me3_peaks$X), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
H3K27me3_peaks <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=H3K27me3_peaks$cluster)
H3K27me3_peaks$start <- as.numeric(H3K27me3_peaks$start)
H3K27me3_peaks$end <- as.numeric(H3K27me3_peaks$end)
H3K27me3_peaks <- as.data.table(H3K27me3_peaks)
H3K27me3_peaks$H3K27me3_label <- paste0(H3K27me3_peaks$chr,":",H3K27me3_peaks$start,"-",H3K27me3_peaks$end)
setDT(H3K27me3_peaks)
setkey(H3K27me3_peaks,chr,start,end)
overlaps <- foverlaps(PMD_HMD_region, H3K27me3_peaks, type = "any", nomatch = 0L)  
overlaps <- as.data.frame(overlaps)

to_plot <- overlaps %>%
  group_by(cluster, V5) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(cluster) %>%
  mutate(total_count = sum(count),
         proportion = count / total_count * 100)
to_plot <- as.data.frame(to_plot)
to_plot$cluster <- paste0("kmeans",to_plot$cluster)
ggplot(to_plot, aes(x = cluster, y = proportion, fill = V5)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle("H3K27me3 peaks")
