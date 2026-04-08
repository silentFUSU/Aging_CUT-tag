rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
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
tissue <- "ovary"
DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01.txt"),header = T)
DMR <- as.data.table(DMR)
setDT(DMR)
setkey(DMR,chr,start,end)
hyper <- DMR[which(DMR$areaStat > 0),c("chr","start","end","meanMethy1","meanMethy2","nCG","areaStat")]
hypo <- DMR[which(DMR$areaStat < 0),c("chr","start","end","meanMethy1","meanMethy2","nCG","areaStat")]

DAR <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
DAR <- as.data.table(DAR)
setDT(DAR)
setkey(DAR,Chr,Start,End)
DAR_Up <- DAR[which(DAR$Significant=="Up"),c("Chr","Start","End")]
DAR_Down <- DAR[which(DAR$Significant=="Down"),c("Chr","Start","End")]


hyper_Down <- foverlaps(hyper, DAR_Down, type = "any", nomatch = 0L)
hyper_Down <- hyper_Down[,c("chr","start","end","meanMethy1","meanMethy2","nCG","areaStat")]
write.table(hyper_Down,paste0("data/samples/WGBS/",tissue,"/DSS_table/bed/",tissue,"_DMR_hyper_DAR_Down.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
hypo_Up <- foverlaps(hypo, DAR_Up, type = "any", nomatch = 0L)
hypo_Up <- hypo_Up[,c("chr","start","end","meanMethy1","meanMethy2","nCG","areaStat")]
write.table(hypo_Up,paste0("data/samples/WGBS/",tissue,"/DSS_table/bed/",tissue,"_DMR_hypo_DAR_Up.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)

