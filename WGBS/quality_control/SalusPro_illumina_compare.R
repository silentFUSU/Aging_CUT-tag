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
    }
  }
  return(tissue_label)
}
SalusPro_path <- "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20241028_SalusPro_WGBS/bed/"
illumina_path <- "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/bonemarrow/bdg/"

samples <- c("DYQ027","DYQ028","DYQ029","DYQ030")
for(sample in samples){
  illumina <- fread(paste0(illumina_path,sample,"_CpG.bdg"),sep = "\t")
  illumina[, percent := V4/V5*100]
  
  saluspro <- fread(paste0(SalusPro_path,"WGBS_bone_marrow_20241028_",sample,"_CpG.bdg"),sep = "\t")
  saluspro[, percent := V4/V5*100]
  
  illumina_sub <- illumina[V5 > 15, .(V1, V2, V3, percent)]  
  saluspro_sub <- saluspro[V5 > 15, .(V1, V2, V3, percent)]  
  
  merged_data <- merge(  
    illumina_sub,   
    saluspro_sub,   
    by = c("V1", "V2", "V3"),   
    suffixes = c(".illumina", ".saluspro")  
  )  
  merged_data <- as.data.frame(merged_data)
  par(cex.lab = 1.5, cex.axis = 1.2)  
  smoothScatter(merged_data[,4] ~ merged_data[,5],xlab = "SalusPro",ylab = "illumina",main = sample,xlim = c(0,100),ylim = c(0,100))
  abline(a = 0, b = 1, col = "red", lty = 2)  
  mtext(paste0(colnames(merged_data)[4]," = ", round(mean(merged_data[,4]),2),"%"), side = 3, line = -2, adj = 0.05,  cex = 1.2)  
  mtext(paste0(colnames(merged_data)[5]," = ", round(mean(merged_data[,5]),2),"%"), side = 3, line = -3.5, adj = 0.05,  cex = 1.2)  
  mtext(paste0("r = ", round(cor(merged_data[,2],merged_data[,3]),2)), side = 1, line = -1.5, adj = 0.9,  cex = 1.2)  
  
}
