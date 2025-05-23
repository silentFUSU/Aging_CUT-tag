rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
peak_list <- data.frame()
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K9me3"
window_size <- "5000"
gap_size <- "10000"
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
peak_pool <- data.frame()
for(tissue in tissues){
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100.bed"))
  peaks$length <- peaks$V3 - peaks$V2 +1
  peaks$tissue <- tissue_label_change(tissue)
  peak_pool <- rbind(peak_pool,peaks)
}
peak_pool$label <- paste(peak_pool$V1,peak_pool$V2,peak_pool$V3,sep = "-")
peak_pool <- peak_pool %>%  
  distinct(label, .keep_all = TRUE)   
# ggplot(peak_pool, aes(x = log10(length))) +  
#   geom_histogram(binwidth = 0.1, color = "black", fill = "skyblue") +  
#   labs(title = "Distribution of Length", x = "log10(Length)", y = "Frequency") +  
#   theme_minimal() #确认peak过滤指标是50kb
peak_pool <- peak_pool[which(peak_pool$length >=50000),]
# write.table(peak_pool[,1:3], paste0("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_peaks_pool.bed"),sep = "\t",col.names = F,row.names = F,append = F,quote = F)
peak_pool$label <- paste0("peak",c(1:nrow(peak_pool)))
peak_pool$V2 <- peak_pool$V2 + 1
peaks_set <- data.frame()
while(nrow(peak_pool) != 0){
  peak_pool <- peak_pool[order(-peak_pool$length, peak_pool$V1, peak_pool$V2),] # 寻找最大length的peak
  max_peak <- peak_pool[1,c(1:4)]
  peak_pool <- peak_pool[-1,]
  check=0
  while(check != 1){
    max_peak <- as.data.table(max_peak)
    setDT(max_peak)
    setkey(max_peak,V1,V2,V3)
    t_peak_pool <- as.data.table(peak_pool)
    setDT(t_peak_pool)
    setkey(t_peak_pool,V1,V2,V3)
    overlaps <- foverlaps(t_peak_pool, max_peak, type = "any", nomatch = 0L)  
    overlaps <- as.data.frame(overlaps)
    if(nrow(overlaps)==0){
      check=1
    }else{
      overlaps$overlap_length <- pmax(0, pmin(overlaps$V3, overlaps$i.V3) - pmax(overlaps$V2, overlaps$i.V2))  
      overlaps <- overlaps[which(overlaps$overlap_length >= overlaps$i.length/2 | overlaps$overlap_length >= overlaps$length/2),]
      if(nrow(overlaps)==0){
        check=1
      }else{
        peak_pool <- peak_pool[-which(peak_pool$label %in% overlaps$label),]
        max_peak <- data.frame(V1=max_peak$V1,V2=min(overlaps$i.V2, max_peak$V2),V3=max(overlaps$i.V3, max_peak$V3))
        max_peak$length <- max_peak$V3 - max_peak$V2 + 1
      }
    }
  }
  peaks_set <- rbind(peaks_set,max_peak[,c(1:4)])
  print(nrow(peaks_set))
  print(nrow(peak_pool))
}
peaks_set <- peaks_set[order(peaks_set$V1,peaks_set$V2,peaks_set$V3),]
write.table(peaks_set[,1:3], paste0("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed"),sep = "\t",col.names = F,row.names = F,append = F,quote = F)

ggplot(peaks_set, aes(x = log10(length))) +
  geom_histogram(binwidth = 0.1, color = "black", fill = "skyblue") +
  labs(title = "Distribution of Length", x = "log10(Length)", y = "Frequency") +
  theme_minimal() #最终peak set的片段长度分布

