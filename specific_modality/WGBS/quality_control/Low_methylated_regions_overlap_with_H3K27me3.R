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
library(dplyr)
library(GenomeInfoDb)
library("GenomicRanges")
library(genomation)
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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/10kb_bins_all_depth.csv"))
  df <- df[which(df$total_V5 > 15),]
  average_percent_by_label <- df %>%
    group_by(label) %>%
    summarise(mean_percent = mean(percent))
  average_percent_by_label <- average_percent_by_label[order(average_percent_by_label$mean_percent),]
  # average_percent_by_label <- average_percent_by_label[c(1:1000),]
  average_percent_by_label <- average_percent_by_label[average_percent_by_label$mean_percent <0.2,]
  bin <- as.data.table(read.table("~/ref_data/mm10_10kb_bins.bed"))
  setDT(bin)
  setkey(bin,V1,V2,V3)
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peak <- peak[,c(1:3)]
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  overlaps <- foverlaps(peak, bin, type = "any", nomatch = 0L)  
  average_percent_by_label$condition <- "out"
  average_percent_by_label$condition[which(average_percent_by_label$label%in% overlaps$V4)] <- "in peak"
  t_tissue_summary <- as.data.frame(table(average_percent_by_label$condition))
  t_tissue_summary$percent <- t_tissue_summary$Freq/sum(t_tissue_summary$Freq)*100
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}
to_plot <- tissue_summary
tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
to_plot$tissue <- factor(to_plot$tissue,levels=tissue_order)
ggplot(to_plot, aes(x = tissue, y =percent, fill = Var1)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  # scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("percent")+
  ggtitle(NULL) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")

