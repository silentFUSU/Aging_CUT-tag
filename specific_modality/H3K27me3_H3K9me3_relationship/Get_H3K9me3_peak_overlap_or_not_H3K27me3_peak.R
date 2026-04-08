rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE) 
tissue <- args[1]  
H3K9me3_peak <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_young_merge-W1000-G3000-E100_compress.bed"))
H3K9me3_peak$V2 <- H3K9me3_peak$V2+1

overlap_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/H3K9me3_H3K27me3_peak_overlap_region.bed"))
overlap_region$V2 <- overlap_region$V2+1

H3K9me3_peak$length <- H3K9me3_peak$V3 - H3K9me3_peak$V2
overlap_region$length <- overlap_region$V3 - overlap_region$V2

H3K9me3_peak <- as.data.table(H3K9me3_peak)
overlap_region <- as.data.table(overlap_region)

setDT(H3K9me3_peak)
setDT(overlap_region)
setkey(H3K9me3_peak,V1,V2,V3)
setkey(overlap_region,V1,V2,V3)

overlap <- foverlaps(overlap_region, H3K9me3_peak, type = "any", nomatch = 0L) 
result <- overlap %>%  
  group_by(V4) %>%  
  summarise(overlap_length = sum(i.length))  

result <- merge(result,H3K9me3_peak,by="V4")
result$percent <- result$overlap_length/result$length * 100

H3K9me3_peak <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_young_merge-W1000-G3000-E100_compress.bed"))
overlap_peaks <- result$V4[which(result$percent > 20)]
not_overlap_peaks <- H3K9me3_peak[-which(H3K9me3_peak$V4 %in% overlap_peaks),]
overlap_peaks <- H3K9me3_peak[which(H3K9me3_peak$V4 %in% overlap_peaks),]
write.table(not_overlap_peaks,paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/H3K9me3_not_overlap_with_H3K27me3_peaks.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
write.table(overlap_peaks,paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/H3K9me3_overlap_with_H3K27me3_peaks.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)

