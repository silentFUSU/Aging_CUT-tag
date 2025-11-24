rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggsignif)
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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

peaks <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed")
peaks <- peaks[which((peaks$V3-peaks$V2 + 1) > 200000),]
peaks$label <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)

kmeans <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
colnames(kmeans)[1] <- "label"

peaks <- merge(peaks, kmeans, by="label",all=T)
peaks$cluster[is.na(peaks$cluster)] <- "Stable"
peaks$length <- peaks$V3-peaks$V2+1
peaks$cluster[which(peaks$cluster %in% c("1","2","3","4"))] <- paste0("kmeans",peaks$cluster[which(peaks$cluster %in% c("1","2","3","4"))])

peaks <- peaks[,c("label","cluster")]

random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
random_regions1$label <- paste0(random_regions1$V1,":",random_regions1$V2,"-",random_regions1$V3)
peaks <- rbind(peaks,random_regions1[,c("label","cluster")])

random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
random_regions2$label <- paste0(random_regions2$V1,":",random_regions2$V2,"-",random_regions2$V3)
peaks <- rbind(peaks,random_regions2[,c("label","cluster")])
# peaks <- peaks %>%
#   mutate(cluster = ifelse(grepl("chrY", label) & cluster == "kmeans1", "kmeans1-chrY", cluster))


antibody <- "H3K27me3"
window_size=5000
gap_size=10000
summary <- data.frame()
age="3m"
for(tissue in tissues){
  if(antibody == "H3K9me3"){
    df <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion.counts"),skip=1)
    df2 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY.counts"),skip=1)
    df3 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY_out_recursion_peaks.counts"),skip=1)
    df <- rbind(df,df2)
    df <- rbind(df,df3)
    df_summary <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion.counts.summary"),header = T)
  }else{
    df <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts"),skip=1)
    df2 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY.counts"),skip=1)
    df3 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY_out_recursion_peaks.counts"),skip=1)
    df <- rbind(df,df2)
    df <- rbind(df,df3)
    df_summary <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)
  }
  
  if(tissue %in% c("mammarygland","ovary","uterus")){
    df <- df[which(df$Chr %in% paste0("chr",c(1:19,"X"))),]
  }else{
    df[which(df$Chr=="chrY"),7:ncol(df)] <- 2*  df[which(df$Chr=="chrY"),7:ncol(df)]
  }
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
  colnames(df_summary) <- gsub(pattern,"\\1",colnames(df_summary))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$age==age & search_table$tissue==tissue & search_table$antibody==antibody),]
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
  df <- merge(df,peaks[,c("label","cluster")],by="label")
  df$tissue <- tissue_label_change(tissue)
  # df$cluster[which(df$Chr=="chrY" & df$cluster=="kmeans1")] <- "chrY_kmeans1"
  summary <- rbind(summary,df)
}
# write.csv(summary,paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_",antibody,"_signal_detail.csv"))

tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary","Colon","Stomach",
                  "Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
summary$tissue <- factor(summary$tissue,levels=tissue_order)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
summary$cluster <- factor(summary$cluster,levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
if(age=="3m"){
  age_label <- "Young"
}else{
  age_label <- "Old"
}
ggplot(summary, aes(x = cluster, y =mean_RPKM, fill = tissue)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  scale_fill_manual(values = color) +
  labs(  
    title = paste0("Distribution of ",age_label," samples ",antibody," RPKM by Cluster"),  
    x = "Cluster",  
    # y = "log2(mean_RPKM)"  
    y="mean_RPKM"
  ) +  
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45,vjust = 1,hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 

color <- setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(paste0("kmeans",1:4),"Stable"))
ggplot(summary, aes(x = cluster, y =mean_RPKM, fill = cluster)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  scale_fill_manual(values = color) +
  labs(  
    title = paste0("Distribution of Young samples ",antibody," RPKM by Cluster"),  
    x = "Cluster",  
    y="mean RPKM"
  ) + 
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 

result <- summary %>%
  group_by(cluster) %>%
  summarise(median_RPKM = median(mean_RPKM, na.rm = TRUE))
write.csv(result,paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_",antibody,"_signal.csv"))



#### two age boxplot
antibody <- "H3K9me3"
window_size=5000
gap_size=10000
summary <- data.frame()

ages=c("3m","24m")
for(age in ages){
  for(tissue in tissues){
    if(antibody == "H3K9me3"){
      df <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion.counts"),skip=1)
      df2 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY.counts"),skip=1)
      df3 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY_out_recursion_peaks.counts"),skip=1)
      df <- rbind(df,df2)
      df <- rbind(df,df3)
      df_summary <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion.counts.summary"),header = T)
    }else{
      df <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts"),skip=1)
      df2 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY.counts"),skip=1)
      df3 <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks_random_200kb_region_rmchrY_out_recursion_peaks.counts"),skip=1)
      df <- rbind(df,df2)
      df <- rbind(df,df3)
      df_summary <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)
    }
    
    if(tissue %in% c("mammarygland","ovary","uterus")){
      df <- df[which(df$Chr %in% paste0("chr",c(1:19,"X"))),]
    }else{
      df[which(df$Chr=="chrY"),7:ncol(df)] <- 2*  df[which(df$Chr=="chrY"),7:ncol(df)]
    }
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
    colnames(df_summary) <- gsub(pattern,"\\1",colnames(df_summary))
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    search_table <- search_table[which(search_table$age==age & search_table$tissue==tissue & search_table$antibody==antibody),]
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
    df <- merge(df,peaks[,c("label","cluster")],by="label")
    df$tissue <- tissue_label_change(tissue)
    df$age <- age
    # df$cluster[which(df$Chr=="chrY" & df$cluster=="kmeans1")] <- "chrY_kmeans1"
    summary <- rbind(summary,df)
  }
}
result <- summary %>%
  group_by(label,age,cluster) %>%
  summarise(median_RPKM = median(mean_RPKM, na.rm = TRUE))

result$cluster <- factor(result$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
result$age <- factor(result$age, levels = c("3m","24m"))
color <- setNames(c("#f39b7f","#4dbbd5"),c("3m","24m"))
p <- ggplot(result, aes(x = cluster, y =median_RPKM, fill = age)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  scale_fill_manual(values = color) +
  labs(  
    title = paste0("Distribution of ",antibody," RPKM by Cluster"),  
    x = "Cluster",  
    y="RPKM"
  ) + 
  theme_bw() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 
ggsave(paste0("result/Sup_figures/H3K9me3_kmeans_",antibody,"_RPKM_per_age.pdf"),p,height = 4,width = 8)
