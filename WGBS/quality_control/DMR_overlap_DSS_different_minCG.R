rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)  
library(DSS)
library(patchwork)
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

tissues <- c("liver","lung","kidney","ileum","Hip","mammarygland","skin","bonemarrow",
             "jejunum","colon","ovary","CB","BAT","thymus","testis","stomach","heart",
             "muscle","bladder","aorta","tongue","spleen","pancreas","brain",
             "cecum","uterus","iWAT")
DMR_summary_minCG3 <- data.frame(condition = as.character(),
                          count = as.numeric(),
                          percent = as.numeric(),
                          tissue = as.character())
DMR_summary_minCG5 <- data.frame(condition = as.character(),
                              count = as.numeric(),
                              percent = as.numeric(),
                              tissue = as.character())

tissue <- "lung"
overlap_summary <- data.frame()
for(tissue in tissues){
  DMR_minCG3 <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01.txt"),header = T)
  DMR_minCG3 <- DMR_minCG3[,c(1:3,9)]
  DMR_minCG3$label_minCG3 <- paste0(DMR_minCG3$chr,":",DMR_minCG3$start,"-",DMR_minCG3$end)
  DMR_minCG3_increase <- as.data.table(DMR_minCG3[which(DMR_minCG3$areaStat > 0),])
  setDT(DMR_minCG3_increase) 
  setkey(DMR_minCG3_increase,chr,start,end)
  DMR_minCG3_decrease <- as.data.table(DMR_minCG3[which(DMR_minCG3$areaStat < 0),])
  setDT(DMR_minCG3_decrease)
  setkey(DMR_minCG3_decrease,chr,start,end)
  
  
  DMR_minCG5 <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01_minCG5.txt"),header = T)
  DMR_minCG5 <- DMR_minCG5[,c(1:3,9)]
  DMR_minCG5$label_minCG5 <- paste0(DMR_minCG5$chr,":",DMR_minCG5$start,"-",DMR_minCG5$end)
  DMR_minCG5_increase <- as.data.table(DMR_minCG5[which(DMR_minCG5$areaStat > 0),])
  setDT(DMR_minCG5_increase) 
  setkey(DMR_minCG5_increase,chr,start,end)
  DMR_minCG5_decrease <- as.data.table(DMR_minCG5[which(DMR_minCG5$areaStat < 0),])
  setDT(DMR_minCG5_decrease)
  setkey(DMR_minCG5_decrease,chr,start,end)
  
  DMR_increase_overlap <- foverlaps(DMR_minCG5_increase, DMR_minCG3_increase, type = "any", nomatch = 0L)  
  DMR_increase_overlap <-unique(DMR_increase_overlap$label_minCG3)
  
  DMR_decrease_overlap <- foverlaps(DMR_minCG5_decrease,DMR_minCG3_decrease, type = "any", nomatch = 0L)
  DMR_decrease_overlap <- unique(DMR_decrease_overlap$label_minCG3)
  
  t_overlap_summary <- data.frame(condition=c("increase","decrease"),
                                  DSS_count=c(nrow(DMR_minCG3_increase),nrow(DMR_minCG3_decrease)),
                                  overlap_count=c(length(DMR_increase_overlap),length(DMR_decrease_overlap)),
                                  tissue=c(tissue_label_change(tissue),tissue_label_change(tissue)))
  
  overlap_summary <- rbind(overlap_summary,t_overlap_summary)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(overlap_summary$tissue)))
overlap_summary$DSS_count[which(overlap_summary$condition=="decrease")] <- -overlap_summary$DSS_count[which(overlap_summary$condition=="decrease")]
overlap_summary$overlap_count[which(overlap_summary$condition=="decrease")] <- -overlap_summary$overlap_count[which(overlap_summary$condition=="decrease")]

overlap_summary$percent <- abs(overlap_summary$overlap_count) / abs(overlap_summary$DSS_count) * 100
ggplot(overlap_summary, aes(x = tissue, y = DSS_count, fill = tissue,color = condition)) +  
  geom_bar(stat = "identity",aes(alpha = ifelse(condition == "Decrease", 0.8, 1))) + 
  geom_bar(data = overlap_summary,aes(x = tissue, y = overlap_count), stat = "identity", fill = "black", alpha = 0.5) + 
  geom_text(aes(label = paste0(round(percent, 1), "%"), y = DSS_count + 100), size = 3, hjust = 0.5,color="black") +
  theme_minimal() +  
  xlab(NULL)+
  ylab("Count")+
  scale_y_continuous(labels = abs) +  
  scale_fill_manual(values=color) +
  theme_bw()+
  theme(  
    axis.title.x = element_text(size = 14),     
    axis.title.y = element_text(size = 14),    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),   
    axis.text.y = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold"),
  ) +  
  scale_color_manual(values = c("Increase" = "white","Decrease" = "black")) 


