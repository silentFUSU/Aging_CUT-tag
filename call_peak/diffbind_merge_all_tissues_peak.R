rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(DiffBind)
tissue_peak_information <- data.frame(SampleID=as.character(),
                                      tissue=as.character(),
                                      Peaks=as.character(),
                                      PeakCaller=as.character())
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K9me3"
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
for(tissue in tissues){
  t_df <- data.frame(SampleID=tissue_label_change(tissue),
                     tissue=tissue_label_change(tissue),
                     Peaks=paste0("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W1000-G3000-E100.bed"),
                     PeakCaller="bed")
  tissue_peak_information <- rbind(tissue_peak_information,t_df)
}
dba_obj <- dba(sampleSheet=tissue_peak_information)
dba.plotHeatmap(dba_obj)
