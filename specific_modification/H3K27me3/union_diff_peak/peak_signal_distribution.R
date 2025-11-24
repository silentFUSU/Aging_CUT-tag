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

peak_num <- read.csv("data/samples/all/H3K27me3/peaks_merged/merged_peaks_tissue_num.csv",row.names = 1)
summary <- data.frame()
for(tissue in tissues){
  df <- read.delim(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_0_tissues.counts"),skip=1)
  df_summary <- read.delim(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_0_tissues.counts.summary"),header = T)
  if(tissue %in% c("mammarygland","ovary","uterus")){
    df <- df[which(df$Chr %in% paste0("chr",c(1:19,"X"))),]
  }else{
    df[which(df$Chr=="chrY"),7:ncol(df)] <- 2*  df[which(df$Chr=="chrY"),7:ncol(df)]
  }
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
  colnames(df_summary) <- gsub(pattern,"\\1",colnames(df_summary))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$age=="3m" & search_table$tissue==tissue & search_table$antibody==antibody),]
  samples <- search_table$sample_name
  df <- data.frame(df[,c(1:6)],df[,samples])
  df_summary <- df_summary[-2,samples]
  colnames(df)[1] <- "label"
  length_kb <- df$Length / 1000  
  
  total_reads <- colSums(df_summary)
  if(! tissue %in% c("mammarygland","ovary","uterus")){
    total_reads <- total_reads + colSums(df[which(df$Chr=="chrY"),7:ncol(df)])/2
  }
  total_reads_million <- total_reads / 1e6  
  for (i in c(7:ncol(df))) {  
    df[[i]] <- (df[[i]] / (length_kb * total_reads_million[i - 6]))  
  }  
  df$mean_RPKM <- rowMeans(df[,c(7:ncol(df))])
  df <- df[,c("label","Chr","mean_RPKM")]
  df <- merge(df,peak_num[,c("label","num")],by="label")
  df$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,df)
  }
to_plot <- summary %>%
  group_by(tissue, num) %>%
  summarise(average_mean_RPKM = mean(mean_RPKM, na.rm = TRUE))

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$tissue)))
to_plot <- as.data.frame(to_plot)
to_plot$num <- as.character(to_plot$num)
to_plot$num <- factor(to_plot$num,levels = c(1:27))
ggplot(to_plot, aes(x = num, y =average_mean_RPKM,fill=num)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  # scale_fill_manual(values = color) +
  labs(  
    x = "tissue num",  
    # y = "log2(mean_RPKM)"  
    y="mean_RPKM"
  ) +  
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 
