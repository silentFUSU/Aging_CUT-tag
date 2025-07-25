rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
antibody <- "H3K9me3"
peaks_quality_control <- function(antibody){
  length_summary <- data.frame()
  broad_peaks_length_summary <- data.frame()
  count_summary <- data.frame()
  broad_peaks_count_summary <- data.frame()
  tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
               "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
  for(tissue in tissues){
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      window_size="5000"
      gap_size="10000"
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"))
    }
    genome <- read.table("~/ref_data/for_normal_mapping/mm10/mm10.normal.chrom.sizes")
    if(tissue %in% c("ovary","mammarygland","uterus")){
      peaks <- peaks[-which(peaks$V1 == "chrY"),]
      genome_size <- sum(genome$V2[which(genome$V1 %in% paste0("chr",c(1:19,"X")))])
    }else{
      genome_size <- sum(genome$V2)
    }
    peaks$length <- peaks$V3 - peaks$V2 +1
    t_count_summary <- data.frame(tissue=tissue_label_change(tissue),count=nrow(peaks))
    t_broad_peaks_count_summary <- data.frame(tissue=tissue_label_change(tissue),count=nrow(peaks[which(peaks$length > 200000),]))
    t_length_summary <- data.frame(tissue=tissue_label_change(tissue),coverage=sum(peaks$length)/genome_size*100)
    t_broad_peaks_length_summary <- data.frame(tissue=tissue_label_change(tissue),coverage=sum(peaks$length[which(peaks$length>200000)])/genome_size*100)
    count_summary <- rbind(count_summary,t_count_summary)
    broad_peaks_count_summary <- rbind(broad_peaks_count_summary,t_broad_peaks_count_summary)
    length_summary <- rbind(length_summary,t_length_summary)
    broad_peaks_length_summary <- rbind(broad_peaks_length_summary,t_broad_peaks_length_summary)
  }
  count_to_plot <- merge(count_summary,broad_peaks_count_summary,by="tissue")
  colnames(count_to_plot) <- c("tissue","count","broad_peaks_count")
  count_to_plot$label <- paste0(count_to_plot$broad_peaks_count,"/",count_to_plot$count)
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- color$V1
  color <- setNames(color,sort(unique(count_to_plot$tissue)))
  ggplot(count_to_plot,mapping = aes(x=tissue,y=count,fill = tissue))+
    geom_bar(stat = "identity", position = position_dodge2())+
    theme_bw()+ylab("")+
    geom_bar(aes(y = broad_peaks_count), stat = "identity", fill = "black", alpha = 0.5) + 
    ylab(antibody)+
    xlab(NULL)+
    ggtitle(paste0(antibody," peaks (broad/all)"))+
    theme(  
      text = element_text(size = 10),  
      axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14), 
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1) 
    ) +
    scale_fill_manual(values = color) +
    theme(legend.position = "none") + 
    geom_text(aes(label = label), position = position_dodge2(width = 0.9), size = 3) 
  
  length_to_plot <- merge(length_summary,broad_peaks_length_summary,by="tissue")
  colnames(length_to_plot) <- c("tissue","coverage","broad_peaks_coverage")
  length_to_plot$label <- paste0(round(length_to_plot$broad_peaks_coverage,2),"/",round(length_to_plot$coverage,2))
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- color$V1
  color <- setNames(color,sort(unique(count_to_plot$tissue)))
  
  ggplot(length_to_plot,mapping = aes(x=tissue,y=coverage,fill = tissue))+
    geom_bar(stat = "identity", position = position_dodge2())+
    theme_bw()+ylab("")+
    geom_bar(aes(y = broad_peaks_coverage), stat = "identity", fill = "black", alpha = 0.5) + 
    ylab(antibody)+
    xlab(NULL)+
    ggtitle(paste0(antibody," peaks coverage (broad/all)"))+
    theme(  
      text = element_text(size = 10),  
      axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 14), 
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1) 
    ) +
    scale_fill_manual(values = color) +
    theme(legend.position = "none") + 
    geom_text(aes(label = label), position = position_dodge2(width = 0.9), size = 3) 
}
