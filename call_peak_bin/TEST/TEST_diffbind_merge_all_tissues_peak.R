rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(DiffBind)
options(scipen = 999)
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
window_size <- "5000"
gap_size <- "10000"
for(tissue in tissues){
  t_df <- data.frame(SampleID=tissue_label_change(tissue),
                     tissue=tissue_label_change(tissue),
                     Peaks=paste0("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"),
                     PeakCaller="bed")
  tissue_peak_information <- rbind(tissue_peak_information,t_df)
}
dba_obj <- dba(sampleSheet=tissue_peak_information)
dba.plotHeatmap(dba_obj)
peaks <- as.data.frame(dba_obj$merged)
peaks$CHR <- paste0("chr",peaks$CHR)
peaks$CHR[which(peaks$CHR=="chr20")] <- "chrX"
peaks$CHR[which(peaks$CHR=="chr21")] <- "chrY"
peaks <- peaks[,c("CHR","START","END")]
peaks$peaks <- paste0("peak",c(1:nrow(peaks)))
write.table(peaks,paste0("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W",window_size,"-G",gap_size,"-E100_diffbind.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
peaks$Length <- peaks$END-peaks$START+1
