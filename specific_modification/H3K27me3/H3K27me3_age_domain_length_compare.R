rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
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
tissues <-  c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
              "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

summary <- data.frame()
for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/edd_peaks_fdr05.bed"))
  domain$length <- domain$V3-domain$V2+1
  domain$tissue <- tissue_label_change(tissue)
  domain$condition <- "H3K27me3 increase"
  domain <- domain[,c("length","tissue","condition")]
  summary <- rbind(summary,domain)  
  }

for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/young_edd_peaks_fdr05.bed"))
  domain$length <- domain$V3-domain$V2+1
  domain$tissue <- tissue_label_change(tissue)
  domain$condition <- "H3K27me3 decrease"
  domain <- domain[,c("length","tissue","condition")]
  summary <- rbind(summary,domain)  
}

for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K9me3/peaks/edd/edd_peaks_fdr05.bed"))
  domain$length <- domain$V3-domain$V2+1
  domain$tissue <- tissue_label_change(tissue)
  domain$condition <- "H3K9me3 increase"
  domain <- domain[,c("length","tissue","condition")]
  summary <- rbind(summary,domain)  
}

for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K9me3/peaks/edd/young_edd_peaks_fdr05.bed"))
  domain$length <- domain$V3-domain$V2+1
  domain$tissue <- tissue_label_change(tissue)
  domain$condition <- "H3K9me3 decrease"
  domain <- domain[,c("length","tissue","condition")]
  summary <- rbind(summary,domain)  
}

for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K36me3/peaks/edd/edd_peaks_fdr05.bed"))
  domain$length <- domain$V3-domain$V2+1
  domain$tissue <- tissue_label_change(tissue)
  domain$condition <- "H3K36me3 increase"
  domain <- domain[,c("length","tissue","condition")]
  summary <- rbind(summary,domain)  
}

for(tissue in tissues){
  domain <- read.table(paste0("data/samples/",tissue,"/H3K36me3/peaks/edd/young_edd_peaks_fdr05.bed"))
  domain$length <- domain$V3-domain$V2+1
  domain$tissue <- tissue_label_change(tissue)
  domain$condition <- "H3K36me3 decrease"
  domain <- domain[,c("length","tissue","condition")]
  summary <- rbind(summary,domain)  
}

color_tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                        "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum"))
color_tissues <- sapply(color_tissues, tissue_label_change)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(color_tissues))
ggplot(summary,aes(x=condition,y=length,fill=tissue))+    
  geom_boxplot()+
  scale_fill_manual(values = color)+
  ggtitle(paste0("edd domain length"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+labs(fill = "", color = "") +ylab("Length")

