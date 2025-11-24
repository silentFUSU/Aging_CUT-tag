rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)

search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
domain_list <- data.frame()
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

domain_pool <- data.frame()
for(tissue in tissues){
  domains <- read.table(paste0("data/samples/",tissue,"/",antibody,"/peaks/edd/edd_peaks_fdr05.bed"))
  domains$length <- domains$V3 - domains$V2 +1
  domains$tissue <- tissue_label_change(tissue)
  domain_pool <- rbind(domain_pool,domains)
}

domain_pool$label <- paste(domain_pool$V1,domain_pool$V2,domain_pool$V3,sep = "-")
domain_pool <- domain_pool %>%  
  distinct(label, .keep_all = TRUE)   

ggplot(domain_pool, aes(x = log10(length))) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "skyblue") +
  labs(title = "Distribution of Length", x = "log10(Length)", y = "Frequency") +
  theme_minimal()

domain_pool <- domain_pool[which(domain_pool$length > 300000),]
domain_pool$V1 <- factor(domain_pool$V1,levels=paste0("chr",c(1:19,"X","Y")))
domain_pool$V2 <- as.numeric(domain_pool$V2)
domain_pool$V3 <- as.numeric(domain_pool$V3)
domain_pool <- domain_pool %>%
  arrange(V1, V2, V3)

write.table(domain_pool[,c(1:3)],"data/samples/all/H3K27me3/bed/H3K27me3_edd_domain_pool.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
