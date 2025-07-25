rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
library(ggsignif)
library(GenomicRanges)
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
tissue <- "kidney"
window_size <- "5000"
gap_size <- "10000"
antibody <- "H3K9me3"
gaps <- seq(from=0,to=100000,by=10000)
summary <- data.frame()
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle")
for(tissue in tissues){
  for(gap in gaps){
    histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"))
    histone <- histone[which((histone$V3-histone$V2 +1)>10000),]
    gr <- GRanges(seqnames = histone$V1,
                  ranges = IRanges(start = histone$V2, end = histone$V3))
    merged_gr <- reduce(gr, min.gapwidth = gap)
    merged_histone <- data.frame(
      V1 = as.character(seqnames(merged_gr)),
      V2 = start(merged_gr),
      V3 = end(merged_gr)
    )
    merged_histone <- merged_histone[which((merged_histone$V3-merged_histone$V2 +1)>100000),]
    t_summary <- data.frame(gap=gap,count=nrow(merged_histone),tissue=tissue_label_change(tissue))
    summary <- rbind(summary,t_summary)
  }
}

ggplot(summary, aes(x = gap, y = count, color=tissue)) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE, method = "loess") + 
  geom_line()

gap <- 20000
for(tissue in tissues){
  histone <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"))
  histone <- histone[which((histone$V3-histone$V2 +1)>10000),]
  gr <- GRanges(seqnames = histone$V1,
                ranges = IRanges(start = histone$V2, end = histone$V3))
  merged_gr <- reduce(gr, min.gapwidth = gap)
  merged_histone <- data.frame(
    V1 = as.character(seqnames(merged_gr)),
    V2 = start(merged_gr),
    V3 = end(merged_gr)
  )
  merged_histone <- merged_histone[which((merged_histone$V3-merged_histone$V2 +1)>100000),]
  write.table(merged_histone,paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_",gap,"_merged.bed"),append = F,quote = F,row.names = F,col.names = F,sep = "\t")
}



