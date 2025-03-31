rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(ggplot2)
library(data.table)

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "ATAC"
conditions <- c("Up","Down")
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
tissue_summary <- data.frame()
for(tissue in tissues){
  peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_MJ/all_peaks/",tissue,".fwp.filter.non_overlapping.bed"))
  peaks <- peaks[,c(1:4)]
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  for(condition in conditions){
    diff_peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_MJ/Diff_peaks/Peak_",condition,"_",tissue,".bed"))
    diff_peaks_count <- nrow(diff_peaks)
    diff_peaks_table <- as.data.table(diff_peaks)
    setDT(diff_peaks_table)
    setkey(diff_peaks_table,V1,V2,V3)
    overlaps_young <- foverlaps(diff_peaks_table, peaks, type = "any", nomatch = 0L)  
    overlap_peaks_count <- nrow(overlaps_young)
    t_tissue_summary <- data.frame(tissue=tissue,condition=condition,overlap_peaks=overlap_peaks_count,diff_peaks=diff_peaks_count) 
    tissue_summary <- rbind(tissue_summary,t_tissue_summary)
  }
}
antibodys <- "ATAC"
tissue_summary$peak_label <- paste0(tissue_summary$overlap_peaks,"/",tissue_summary$diff_peaks)
tissue_summary <- tissue_summary %>%  
  mutate(tissue_label = sapply(tissue, tissue_label_change))  
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(tissue_summary$tissue_label)))
for(condition in conditions){
  df <- tissue_summary[which(tissue_summary$condition==condition),]
  # df <- arrange(df, diff_peaks)  
  # df$tissue_label <- factor(df$tissue_label,levels=sort(tissues))
  df$tissue_label <- factor(df$tissue_label,levels=tissues_label)
  p_list <- list()
  p_list[[i]] <-  ggplot(df,mapping = aes(x=diff_peaks,y=tissue_label,fill = tissue_label))+
    geom_bar(stat = "identity", position = position_dodge2())+
    theme_bw()+ylab("")+
    geom_bar(aes(x = overlap_peaks), stat = "identity", fill = "black", alpha = 0.5) + 
    xlab(antibodys[[i]])+
    theme(  
      text = element_text(size = 10),  
      axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
      axis.title.y = element_blank(),  
      axis.ticks.y = element_blank(),  
      legend.position = "none"  
    ) +
    scale_fill_manual(values = color) +
    theme(legend.position = "none") + 
    geom_text(aes(label = peak_label), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) + xlim(0,100000)
  combined_plot <- arrangeGrob(  
    grobs = p_list,  
    ncol = length(p_list),  
    widths = c(1.4,rep(1,(length(p_list)-1))),
    top = textGrob(paste0(condition," peak number"), gp = gpar(fontsize = 15, fontface = "bold"))  
  )  
  grid.draw(combined_plot) 
  # ggsave(paste0("result/all/diff/all_tissues_",condition,"_peak_number_remove_batch_effect.png"), plot = combined_plot, width = 18, height = 6,type="cairo")  
}
