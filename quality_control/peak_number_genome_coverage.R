rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)  

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

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
peaks_summary <- data.frame()
antibody <- "H3K27me3"
window_size = 1000
gap_size = 3000
genome_coverage <- read.table("~/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes")
genome_coverage <- genome_coverage[which(genome_coverage$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
genome_length <- sum(genome_coverage$V2)
for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed")) 
  peaks_number <- nrow(df)
  df$length <- df$V3-df$V2+1
  peaks_coverage <- sum(df$length)
  peaks_coverage <- peaks_coverage / genome_length *100
  large_peaks_coverage <- sum(df$length[which(df$length > 100000)])
  large_peaks_coverage <- large_peaks_coverage / genome_length *100
  t_peaks_summary <- data.frame(tissue = tissue_label_change(tissue),
                                peaks_number=peaks_number,
                                peaks_coverage=peaks_coverage,
                                large_peaks_coverage=large_peaks_coverage)
  peaks_summary <- rbind(peaks_summary,t_peaks_summary)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(peaks_summary$tissue)))
ggplot(peaks_summary,mapping = aes(x=tissue,y=peaks_number,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+
  theme_bw()+ylab("")+
  xlab(NULL)+
  theme(  
    text = element_text(size = 12),  
    axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    legend.position = "none"  
  ) + scale_fill_manual(values = color) +    
  geom_text(aes(label = peaks_number), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3, angle = 90)+
  ylim(0,80000)+
  ggtitle("Peaks Number")

peaks_summary$peaks_label <- paste0(round(peaks_summary$large_peaks_coverage,1),"% / ",round(peaks_summary$peaks_coverage,1),"%")
ggplot(peaks_summary,mapping = aes(x=tissue,y=peaks_coverage,fill = tissue))+
  geom_bar(stat = "identity", position = position_dodge2())+
  geom_bar(aes(y= large_peaks_coverage), stat = "identity", fill = "black", alpha = 0.5) + 
  theme_bw()+ylab("Genome Coverge(%)")+
  xlab(NULL)+
  theme(  
    text = element_text(size = 12),  
    axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    legend.position = "none"  
  ) + scale_fill_manual(values = color) +    
  ylim(0,100)+
  geom_text(aes(label = peaks_label), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3, angle = 90)+
  ggtitle("Peaks Coverage")

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
peaks_summary <- data.frame()
antibody <- "H3K9me3"
window_size = 5000
gap_size = 10000
for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed")) 
  df$length <- df$V3-df$V2+1
  df <- df[,"length",drop=F]
  df$tissue <- tissue_label_change(tissue)
  peaks_summary <- rbind(peaks_summary,df)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(peaks_summary$tissue)))
ggplot(peaks_summary, aes(x = tissue, y = log10(length), fill= tissue)) +  
  geom_violin(adjust = 2.5) +          
  scale_fill_manual(values = color) +    
  geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  
  theme_minimal()+
  ggtitle("Peak Length")+
  theme(  
    text = element_text(size = 12),  
    axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    legend.position = "none"  
  )+
  labs(x = NULL,y = "log10(Peak Length)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
