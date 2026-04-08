rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(DOSE)
library(clusterProfiler)
options(scipen = 0)
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
for(tissue in tissues){
  df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/",tissue,"_DARs.txt"))
  peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
  rownames(peaks) <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
  df <- merge(df,peaks,by="row.names")
  stable <- df[which(df$FDR > 0.9 & abs(df$logFC) < 0.05),c("V1","V2","V3")]
  increase <- df[which(df$FDR < 0.05 & df$logFC > 0),c("V1","V2","V3")]
  decrease <- df[which(df$FDR < 0.05 & df$logFC < 0),c("V1","V2","V3")]
  dir.create("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/stable/",showWarnings = F,recursive = T)
  dir.create("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/up/",showWarnings = F,recursive = T)
  dir.create("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/down/",showWarnings = F,recursive = T)
  write.table(stable,paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/stable/",tissue,".bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(increase,paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/up/",tissue,".bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(decrease,paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/down/",tissue,".bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
}



