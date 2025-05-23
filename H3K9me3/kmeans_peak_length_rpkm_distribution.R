rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
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
peaks <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed")
peaks$label <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)

kmeans <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
colnames(kmeans)[1] <- "label"

peaks <- merge(peaks, kmeans, by="label")
peaks$length <- peaks$V3-peaks$V2+1
peaks$cluster <- paste0("kmeans",peaks$cluster)

ggplot(peaks, aes(x = cluster, y = log10(length), fill = cluster)) +  
  geom_boxplot(alpha = 0.7) +  
  labs(  
    title = "Distribution of Lengths by Cluster",  
    x = "Cluster",  
    y = "log10(Length)"  
  ) +  
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")  
  )  

proportion_data <- peaks %>%  
  group_by(cluster, V1) %>%  
  summarize(count = n()) %>%  
  ungroup() %>%  
  group_by(cluster) %>%  
  mutate(proportion = count / sum(count)) 

proportion_data$V1 <- factor(proportion_data$V1, levels=paste0("chr",c(1:19,"X","Y")))
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,paste0("chr",c(1:19,"X","Y")))
# Create a stacked bar chart  
ggplot(proportion_data, aes(x = cluster, y = proportion, fill = V1)) +  
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = color)+
  labs(  
    title = "Kmeans Cluster Chromosome proportion",  
    x = "Cluster",  
    y = "Proportion"  
  ) +  
  theme_minimal() +  
  scale_y_continuous(labels = scales::percent_format()) +  # Format y-axis as percentage  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  )  




tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
summary <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts"),header = T)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$age=="3m" & search_table$tissue==tissue & search_table$antibody=="H3K9me3"),]
  samples <- search_table$sample_name
  df <- data.frame(df[,c(1:6)],df[,samples])
  colnames(df)[1] <- "label"
  length_kb <- df$Length / 1000  
  total_reads <- colSums(df[, c(7:ncol(df))])
  total_reads_million <- total_reads / 1e6  
  for (i in c(7:ncol(df))) {  
    df[[i]] <- (df[[i]] / (length_kb * total_reads_million[i - 6]))  
  }  
  df$mean_RPKM <- rowMeans(df[,c(7:ncol(df))])
  df <- df[,c("label","mean_RPKM")]
  df <- merge(df,kmeans,by="label")
  df$cluster <- paste0("kmeans",df$cluster)
  summary <- rbind(summary,df)
}
ggplot(summary, aes(x = cluster, y =log2(mean_RPKM), fill = cluster)) +  
  geom_boxplot(alpha = 0.7) +  
  labs(  
    title = "Distribution of Young samples RPKM by Cluster",  
    x = "Cluster",  
    y = "log2(mean_RPKM)"  
  ) +  
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")  
  )  
summary <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts"),header = T)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$age=="3m" & search_table$tissue==tissue & search_table$antibody=="H3K9me3"),]
  samples <- search_table$sample_name
  df <- data.frame(df[,c(1:6)],df[,samples])
  colnames(df)[1] <- "label"
  length_kb <- df$Length / 1000  
  total_reads <- colSums(df[, c(7:ncol(df))])
  total_reads_million <- total_reads / 1e6  
  for (i in c(7:ncol(df))) {  
    df[[i]] <- (df[[i]] / (length_kb * total_reads_million[i - 6]))  
  }  
  df$mean_RPKM <- rowMeans(df[,c(7:ncol(df))])
  df <- df[,c("label","mean_RPKM")]
  colnames(df)[2] <- tissue_label_change(tissue)
  if(nrow(summary)==0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="label")
  }
}
summary <- merge(summary,kmeans,by="label")
rownames(summary) <- summary$label
summary <- summary[,-1]
summary_sorted <- summary[order(summary$cluster),]  
data_for_heatmap <- summary_sorted[, -ncol(summary_sorted)]  

annotation <- summary[,"cluster",drop=F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("white", "red"))(100)  
breaks <- c(seq(0, 10, length.out = 100))  
data_for_heatmap <- data_for_heatmap[,c("Kidney","Muscle","Skin","Liver","Aorta","Testis","Cortex","Tongue","Uterus","Bladder","Stomach","Heart","Hippocampus","Cerebellum","BAT","Lung","Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")]
pheatmap::pheatmap(data_for_heatmap,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, annotation_row = annotation, color = color_palette, clustering_distance_cols="manhattan")

