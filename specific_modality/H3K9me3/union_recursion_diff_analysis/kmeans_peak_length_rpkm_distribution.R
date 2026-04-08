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
peaks <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed")
peaks <- peaks[which((peaks$V3-peaks$V2 + 1) > 200000),]
peaks$label <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)

kmeans <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
colnames(kmeans)[1] <- "label"

peaks <- merge(peaks, kmeans, by="label",all=T)
peaks$cluster[is.na(peaks$cluster)] <- "Stable"
peaks$length <- peaks$V3-peaks$V2+1
peaks$cluster[which(peaks$cluster %in% c("1","2","3","4"))] <- paste0("kmeans",peaks$cluster[which(peaks$cluster %in% c("1","2","3","4"))])
color <- setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(paste0("kmeans",1:4),"Stable"))
p <- ggplot(peaks, aes(x = cluster, y = log10(length), fill = cluster)) +  
  geom_boxplot(alpha = 0.7) +  
  scale_fill_manual(values = color) +
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
ggsave("result/Sup_figures/H3K9me3_kmeans_peak_length.pdf",p,width = 6,height = 8)
peaks$cluster[which(peaks$V1=="chrY" & peaks$cluster=="kmeans1")] <-"kmeans1-chrY"
result <- peaks %>%
  group_by(cluster) %>%
  summarise(median_length = median(length, na.rm = TRUE))
write.csv(result,"data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_length.csv")

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
p<-ggplot(proportion_data, aes(x = cluster, y = proportion, fill= V1)) +  
  geom_bar(stat = "identity", position = "fill",color = "black") + 
  scale_fill_manual(values = color)+
  labs(  
    title = "Kmeans Cluster Chromosome proportion",  
    x = "Cluster",  
    y = "Proportion",
    fill = "Chromosome"
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
# ggsave("result/Sup_figures/H3K9me3_kmeans_chromosome.pdf",p,width = 6,height = 8)
proportion_data$cluster <- factor(proportion_data$cluster,levels=rev(c("kmeans1","kmeans2","kmeans3","kmeans4","Stable")))
p<-ggplot(proportion_data, aes(x = proportion, y = cluster, fill= V1)) +  
  geom_bar(stat = "identity", position = "fill",color = "black") + 
  scale_fill_manual(values = color)+
  labs(  
    title = "Kmeans Cluster Chromosome proportion",  
    x = "Proportion",  
    y = "Cluster",
    fill = "Chromosome"
  ) +  
  theme_bw() +   
  scale_x_continuous(labels = scales::percent_format()) + 
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  )  
p
ggsave("result/Sup_figures/H3K9me3_kmeans_chromosome_ppt.pdf",p,width = 8,height = 6)

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
summary <- data.frame()
antibody <- "H3K27me3"
for(tissue in tissues){
  if(antibody == "H3K9me3"){
    df <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts"),header = T)
  }else{
    df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts"),header = T)
  }
  
  if(tissue %in% c("mammarygland","ovary","uterus")){
    df <- df[-which(df$Chr == "chrY"),]
  }

  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$age=="3m" & search_table$tissue==tissue & search_table$antibody==antibody),]
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
  df <- merge(df,peaks[,c("label","cluster")],by="label")
  summary <- rbind(summary,df)
}
color <- setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(paste0("kmeans",1:4),"Stable"))
p <- ggplot(summary, aes(x = cluster, y =mean_RPKM, fill = cluster)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  scale_fill_manual(values = color) +
  labs(  
    title = paste0("Distribution of Young samples ",antibody," RPKM by Cluster"),  
    x = "Cluster",  
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
  # + ylim(-11,11)
 
# ggsave("result/Sup_figures/H3K9me3_kmeans_H3K9me3_RPKM.pdf",p,width = 6,height = 8)

### H3K27me3/H3K9me3 signal per tissue
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K27me3"
p_list <- list()
for(tissue in tissues){
  # df <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts"),header = T)
  df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts"),header = T)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$age=="3m" & search_table$tissue==tissue & search_table$antibody==antibody),]
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
  df <- merge(df,peaks[,c("label","cluster")],by="label")
  p_list[[tissue]] <- ggplot(df, aes(x = cluster, y = mean_RPKM, fill = cluster)) +  
    geom_boxplot(alpha = 0.7,outliers = F) +  
    scale_fill_manual(values = color) +
    labs(  
      title = tissue_label_change(tissue),  
      x = "Cluster",  
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
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 7)
ggsave("result/figures/H3K27me3_signal_in_H3K9me3_recursion_peaks_per_tissues.png",combined_plot,width = 20,height = 20,type="cairo")


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

