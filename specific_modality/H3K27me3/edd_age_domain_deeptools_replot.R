rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
library(ggsignif)
library(deepToolsDownstream)

tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Frontal Cortex"
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
tissue <- "brain"
antibodys<- c("H3K27me3","H3K9me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac")
for(antibody in antibodys){
  se <- importCount(paste0("result/all/H3K27me3_domain/matrix_merge/",tissue,"_",antibody,"_change_in_H3K27me3_domain.mat.gz"))
  p <- plotProfile(se)+
    theme_bw()+   
    geom_line(linewidth = 1.2)+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))
  ggsave(paste0("result/all/H3K27me3_domain/plot_merge/",tissue,"_",antibody,"_change_in_H3K27me3_domain_replot.pdf"),p,width = 6,height = 4)
}
# "#f7fbff","#08306b"
# "#ffffcc""#800026"
# "#f2f0f7""#3f007d"
# "#fcfbfd","#1f0033"
# 
# 
# '#fcfbfd','#9e9ac8','#2d004b'
