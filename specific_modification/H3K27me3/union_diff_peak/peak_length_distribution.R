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
peak_pool_number_each_tissue <- as.data.frame(table(peak_pool$tissue))
to_plot <- peak_pool_number_each_tissue
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(to_plot$Var1))
ggplot(to_plot, aes(x = Var1, y = Freq,fill=Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = color) +
  labs(
    title = "Barplot of Frequencies",
    x = "Category",
    y = "Frequency"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

peak_pool$label <- paste(peak_pool$V1,peak_pool$V2,peak_pool$V3,sep = "-")
peak_pool <- peak_pool %>%  
  distinct(label, .keep_all = TRUE)   

ggplot(peak_pool, aes(x = log10(length))) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "skyblue") +
  labs(title = "Distribution of Length", x = "log10(Length)", y = "Frequency") +
  theme_minimal()

peak_pool <- peak_pool[which(peak_pool$length <= 400000),]# 100000 -> 300000
peak_pool$V1 <- factor(peak_pool$V1,levels=paste0("chr",c(1:19,"X","Y")))
peak_pool$V2 <- as.numeric(peak_pool$V2)
peak_pool$V3 <- as.numeric(peak_pool$V3)
peak_pool <- peak_pool %>%
  arrange(V1, V2, V3)

write.table(peak_pool[,c(1:3)],"data/samples/all/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100_peak_pool.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

