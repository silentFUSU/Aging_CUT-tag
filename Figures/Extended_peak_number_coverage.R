rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(readr) 
library(ggrepel)
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
antibodys <- c("H3K27ac","H3K4me1","H3K4me3")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
peak_summary <- data.frame()
for(antibody in antibodys){
  for(tissue in tissues){
    df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_young_old_narrowpeak.bed"))    
    df$length <- df$V3-df$V2+1
    t_peak_summary <- data.frame(tissue=tissue,antibody=antibody,count=nrow(df),coverage=sum(df$length))
    peak_summary <- rbind(peak_summary,t_peak_summary)
  }  
}
antibody <- "ATAC"
for(tissue in tissues){
  df <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_macs_young_old_narrowpeak_summits_spm3.bed"))    
  df$length <- df$V3-df$V2+1
  t_peak_summary <- data.frame(tissue=tissue,antibody=antibody,count=nrow(df),coverage=sum(df$length))
  peak_summary <- rbind(peak_summary,t_peak_summary)
}  

to_plot <- peak_summary
to_plot$antibody <- factor(to_plot$antibody,levels=c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
p <- ggplot(to_plot, aes(x = antibody, y = count)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Peak count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  ) + 
  ylim(10000,150000)
ggsave("result/Sup_figures/active_marks_peak_count.pdf",p,height=6,width = 6)

to_plot$coverage <- to_plot$coverage/1000000
p <-ggplot(to_plot, aes(x = antibody, y = coverage)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Peak genome coverage (Mb)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )+ylim(0,120)

ggsave("result/Sup_figures/active_marks_peak_coverage.pdf",p,height=6,width = 6)

antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
peak_summary <- data.frame()
for(antibody in antibodys){
  for(tissue in tissues){
    df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W5000-G10000-E100.bed"))    
    df$length <- df$V3-df$V2+1
    t_peak_summary <- data.frame(tissue=tissue,antibody=antibody,count=nrow(df),coverage=sum(df$length))
    peak_summary <- rbind(peak_summary,t_peak_summary)
  }  
}
to_plot <- peak_summary
to_plot$antibody <- factor(to_plot$antibody,levels=c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
p <-ggplot(to_plot, aes(x = antibody, y = count)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2) +
  # scale_fill_manual(values = color) +
  labs(y = "Peak count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )+ylim(5000,25000)
ggsave("result/Sup_figures/broad_marks_peak_count.pdf",p,height=6,width = 5)

to_plot$coverage <- to_plot$coverage/1000000
p <- ggplot(to_plot, aes(x = antibody, y = coverage)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "Peak genome coverage (Mb)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )+ylim(350,1300)
ggsave("result/Sup_figures/broad_marks_peak_coverage.pdf",p,height=6,width = 5)
