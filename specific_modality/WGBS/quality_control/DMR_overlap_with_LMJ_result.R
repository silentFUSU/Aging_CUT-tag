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
DMR_summary <- data.frame(condition = as.character(),
                          count = as.numeric(),
                          percent = as.numeric(),
                          tissue = as.character())
DMR_summary_LMJ <- data.frame(condition = as.character(),
                          count = as.numeric(),
                          percent = as.numeric(),
                          tissue = as.character())

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01_minCG5.txt"),header = T)
  increase <- nrow(DMR[which(DMR$areaStat > 0),])
  decrease <- nrow(DMR[which(DMR$areaStat < 0),])
  tissue_label <- tissue_label_change(tissue)
  t_DMR_percent <- data.frame(condition = c("Increase","Decrease"),
                              count = c(increase,decrease),
                              percent = c(increase/nrow(DMR)*100, decrease/nrow(DMR)*100),
                              tissue = c(tissue_label,tissue_label)) 
  DMR_summary <- rbind(DMR_summary,t_DMR_percent)
}


for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  DMR_increase <- read.table(paste0("data/samples/WGBS/LMJ_DMR_bed/",tissue,"_hyper.bed"),header = F)
  DMR_decrease <- read.table(paste0("data/samples/WGBS/LMJ_DMR_bed/",tissue,"_hypo.bed"),header = F)
  increase <- nrow(DMR_increase)
  decrease <- nrow(DMR_decrease)
  tissue_label <- tissue_label_change(tissue)
  t_DMR_percent <- data.frame(condition = c("Increase","Decrease"),
                              count = c(increase,decrease),
                              percent = c(increase/nrow(DMR)*100, decrease/nrow(DMR)*100),
                              tissue = c(tissue_label,tissue_label)) 
  DMR_summary_LMJ <- rbind(DMR_summary_LMJ,t_DMR_percent)
}

p1 <- ggplot(DMR_summary, aes(x = tissue, y = count, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ylim(0,500000)+
  ggtitle("DMR")

p2 <- ggplot(DMR_summary_LMJ, aes(x = tissue, y = count, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Counts")+
  ylim(0,500000)+
  ggtitle("DMR_from_LMJ")
p1+p2

tissue <- "lung"
overlap_summary <- data.frame()
for(tissue in tissues){
  DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01_minCG5.txt"),header = T)
  DMR <- DMR[,c(1:3,9)]
  DMR$label <- paste0(DMR$chr,":",DMR$start,"-",DMR$end)
  DMR_increase <- as.data.table(DMR[which(DMR$areaStat > 0),])
  setDT(DMR_increase) 
  setkey(DMR_increase,chr,start,end)
  DMR_decrease <- as.data.table(DMR[which(DMR$areaStat < 0),])
  setDT(DMR_decrease)
  setkey(DMR_decrease,chr,start,end)
  
  DMR_increase_LMJ <- read.table(paste0("data/samples/WGBS/LMJ_DMR_bed/",tissue,"_hyper.bed"),header = F)
  DMR_increase_LMJ <- as.data.table(DMR_increase_LMJ)
  setDT(DMR_increase_LMJ)
  setkey(DMR_increase_LMJ,V1,V2,V3)
  
  DMR_decrease_LMJ <- read.table(paste0("data/samples/WGBS/LMJ_DMR_bed/",tissue,"_hypo.bed"),header = F)
  DMR_decrease_LMJ <- as.data.table(DMR_decrease_LMJ)
  setDT(DMR_decrease_LMJ)
  setkey(DMR_decrease_LMJ,V1,V2,V3)
  
  DMR_increase_overlap <- foverlaps(DMR_increase_LMJ, DMR_increase, type = "any", nomatch = 0L)  
  DMR_increase_overlap <-unique(DMR_increase_overlap$label)
  
  DMR_decrease_overlap <- foverlaps(DMR_decrease_LMJ,DMR_decrease, type = "any", nomatch = 0L)
  DMR_decrease_overlap <- unique(DMR_decrease_overlap$label)
  
  t_overlap_summary <- data.frame(condition=c("increase","decrease"),
                                DSS_count=c(nrow(DMR_increase),nrow(DMR_decrease)),
                                overlap_count=c(length(DMR_increase_overlap),length(DMR_decrease_overlap)),
                                tissue=c(tissue_label_change(tissue),tissue_label_change(tissue)))
  
  overlap_summary <- rbind(overlap_summary,t_overlap_summary)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(overlap_summary$tissue)))
overlap_summary$DSS_count[which(overlap_summary$condition=="decrease")] <- -overlap_summary$DSS_count[which(overlap_summary$condition=="decrease")]
overlap_summary$overlap_count[which(overlap_summary$condition=="decrease")] <- -overlap_summary$overlap_count[which(overlap_summary$condition=="decrease")]

ggplot(overlap_summary, aes(x = tissue, y = DSS_count, fill = tissue,color = condition)) +  
  geom_bar(stat = "identity",aes(alpha = ifelse(condition == "Decrease", 0.8, 1))) + 
  geom_bar(data = overlap_summary,aes(x = tissue, y = overlap_count), stat = "identity", fill = "black", alpha = 0.5) + 
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

