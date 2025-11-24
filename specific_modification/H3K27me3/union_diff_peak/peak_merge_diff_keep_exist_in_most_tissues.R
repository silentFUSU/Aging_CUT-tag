rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)

search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K27me3"

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
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_merge-W5000-G10000-E100.bed"))
  peaks$length <- peaks$V3 - peaks$V2 +1
  peaks$tissue <- tissue_label_change(tissue)
  peak_pool <- rbind(peak_pool,peaks)
}


merge_peak <- read.table("data/samples/all/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100_peak_merged.bed")

merge_peak_final <- data.frame()
progress_bar <- txtProgressBar(min = 0, max = nrow(merge_peak), style = 3)
tissue_num <- 27
tissue_num_summary <- data.frame()
for(i in c(1:nrow(merge_peak))){
  peak <- merge_peak[i,]
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  t_peak_pool <- peak_pool[which(peak_pool$V1 == peak$V1),]
  t_peak_pool <- as.data.table(t_peak_pool)
  setDT(t_peak_pool)
  setkey(t_peak_pool,V1,V2,V3)
  overlaps <- foverlaps(peak, t_peak_pool, type = "any", nomatch = 0L)  
  
  tissue_check <- as.data.frame(table(overlaps$tissue))
  t_tissue_num_summary <- data.frame(label=paste0(peak$V1,":",peak$V2,"-",peak$V3),num=nrow(tissue_check))
  tissue_num_summary <- rbind(tissue_num_summary,t_tissue_num_summary)
  if(nrow(tissue_check)>= tissue_num){
    merge_peak_final <- rbind(merge_peak_final, merge_peak[i,])
  }
  setTxtProgressBar(progress_bar, i)
}
write.csv(tissue_num_summary,"data/samples/all/H3K27me3/peaks_merged/merged_peaks_tissue_num.csv")
to_plot <- as.data.frame(table(tissue_num_summary$num))
ggplot(to_plot, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity",fill="skyblue") +
  theme_minimal() +
  labs(
    title = "Barplot of Frequencies",
    x = "Category",
    y = "Frequency"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label = Freq), vjust = -0.5)

write.table(merge_peak_final,paste0("data/samples/all/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
