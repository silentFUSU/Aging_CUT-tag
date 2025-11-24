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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum"))

summary <- data.frame()
for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/edd_peaks_fdr05.bed"))
  # peaks <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100.bed"))
  peaks <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed")
  domain$overlap <- NA
  
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  for(i in c(1:nrow(domain))){
    t_domain <- domain[i,c(1:3)]
    length <- t_domain$V3-t_domain$V2+1
    t_domain <- as.data.table(t_domain)
    setDT(t_domain)
    setkey(t_domain,V1,V2,V3)
    overlaps <- as.data.frame(foverlaps(t_domain, peaks, type = "any", nomatch = 0L))
    if(nrow(overlaps)>0){
      overlaps$overlap_length <- NA  
      for(j in c(1:nrow(overlaps))){
        if(overlaps[j,"V2"] <= overlaps[j,"i.V2"]){
          left <- as.numeric(overlaps[j,"i.V2"])
        }else{
          left <- as.numeric(overlaps[j,"V2"])
        }
        if(overlaps[j,"V3"] >= overlaps[j,"i.V3"]){
          right <- as.numeric(overlaps[j,"i.V3"])
        }else{
          right <- as.numeric(overlaps[j,"V3"])
        }
        overlaps[j,"overlap_length"] <- right - left + 1
      }
      overlap_length <- sum(overlaps$overlap_length)
      domain[i,"overlap"] <- overlap_length/length*100
    }else{
      domain[i,"overlap"] <- 0
    }
  }
  t_summary <- data.frame(tissue=tissue_label_change(tissue),
                          overlap_peak=sum(domain$overlap>60)/nrow(domain) *100,
                          outside_peak=100-sum(domain$overlap>60)/nrow(domain) *100)
  summary <- rbind(summary,t_summary)
}
summary <- summary[order(summary$overlap_peak),]
to_plot <- reshape2::melt(summary)

color <- setNames(c("#009980","#838B8B"),c("overlap_peak","outside_peak"))
to_plot$tissue <- factor(to_plot$tissue,summary$tissue)
ggplot(to_plot, aes(x = tissue, y = value, fill = variable)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle(paste0("H3K27me3 age-domain overlap with H3K9me3"))



