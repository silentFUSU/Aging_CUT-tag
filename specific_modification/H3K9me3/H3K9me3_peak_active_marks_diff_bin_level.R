rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(data.table)
tissue <- "lung"
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
antibody <- "H3K4me3"
increase_summary <- data.frame()
decrease_summary <- data.frame()
for(tissue in tissues){
  active_mark <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_after_remove_batch_effect.csv"),row.names = 1)
  active_mark$Start <- active_mark$Start + 1
  active_mark <- active_mark[,c("Chr","Start","End","LogFC.old.young","Significant")]
  active_mark <- as.data.table(active_mark)
  setDT(active_mark)
  setkey(active_mark,Chr,Start,End)
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3$Start <- H3K9me3$Start + 1
  H3K9me3_increase <- histone[which(histone$Significant == "Up"),c("Geneid","Chr", "Start", "End")]
  H3K9me3_decrease <- histone[which(histone$Significant == "Down"),c("Geneid","Chr", "Start", "End")]
  
  if(nrow(H3K9me3_increase) > 100){
    H3K9me3_increase <- as.data.table(H3K9me3_increase)
    setDT(H3K9me3_increase)
    setkey(H3K9me3_increase,Chr,Start,End)
    overlaps <- foverlaps(active_mark, H3K9me3_increase, type = "any", nomatch = 0L)  
    
    }
}