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
condition <- "peak"
tissue_num <- 27
if(condition == "domain"){
  regions <- read.table("data/samples/all/H3K27me3/bed/H3K27me3_edd_domain_merged.bed")
  regions$label <- paste0(regions$V1,":",regions$V2,"-",regions$V3)
  kmeans <- read.csv("data/samples/all/H3K27me3/edd_domain_merged/kmeans_annotation_RPKM.csv")
}else{
  regions <- read.table(paste0("data/samples/all/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.bed"))
  regions$label <- paste0(regions$V1,":",regions$V2,"-",regions$V3)
  kmeans <- read.csv(paste0("data/samples/all/H3K27me3/peaks_merged/kmeans_annotation_larger_",tissue_num,"_tissues_RPKM.csv"))
}
colnames(kmeans)[1] <- "label"
regions <- merge(regions, kmeans, by="label")
regions$length <- regions$V3-regions$V2+1
regions$cluster <- paste0("kmeans",regions$cluster)

ggplot(regions, aes(x = cluster, y = length, fill = cluster)) +  
  geom_boxplot(alpha = 0.7,outliers = F) +  
  labs(  
    title = "Distribution of Lengths by Cluster",  
    x = "Cluster",  
    y = "Length"  
  ) +  
  theme_minimal() +  
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")  
  )+
  geom_signif(
    comparisons = list(c("kmeans1", "kmeans2"), c("kmeans1", "kmeans3"),c("kmeans2","kmeans3")),
    textsize = 4,test = "t.test",
    y_position = c(230000, 270000, 250000 ),
    tip_length = c(1/100,1/100,1/100),
    map_signif_level=T
  )

proportion_data <- regions %>%  
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

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
rpkm_summary <- data.frame()
antibody <- "H3K27me3"
for(tissue in tissues){ 
  if(condition=="domain"){
    df <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts.summary"),header = T,row.names = 1)
    summary <- summary[-2,]
  }else{
    df <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.counts.summary"),header = T,row.names = 1)
  }

  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(df)[7:ncol(df)] <-  gsub(pattern, "\\1",colnames(df)[7:ncol(df)])
  colnames(summary) <- gsub(pattern, "\\1",colnames(summary))
  
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$age=="24m" & search_table$tissue==tissue & search_table$antibody==antibody),]
  samples <- search_table$sample_name
  df <- data.frame(df[,c(1:6)],df[,samples])
  summary <- summary[,samples]
  colnames(df)[1] <- "label"
  length_kb <- df$Length / 1000  
  total_reads <- colSums(summary)
  total_reads_million <- total_reads / 1e6  
  for (i in c(7:ncol(df))) {  
    df[[i]] <- (df[[i]] / (length_kb * total_reads_million[i - 6]))  
  }  
  df$mean_RPKM <- rowMeans(df[,c(7:ncol(df))])
  df <- df[,c("label","mean_RPKM")]
  df <- merge(df,kmeans,by="label")
  df$cluster <- paste0("kmeans",df$cluster)
  rpkm_summary <- rbind(rpkm_summary,df)
}
ggplot(rpkm_summary, aes(x = cluster, y =log2(mean_RPKM), fill = cluster)) +  
  geom_boxplot(alpha = 0.7) +  
  labs(  
    title = paste0("Distribution of Old samples ",antibody," RPKM by Cluster"),  
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
  )  + 
  ylim(-10,10)+
  geom_signif(
    comparisons = list(c("kmeans1", "kmeans2"), c("kmeans1", "kmeans3"),c("kmeans2","kmeans3")),
    textsize = 4,test = "t.test",
    y_position = c(7, 9, 8 ),
    tip_length = c(1/100,1/100,1/100),
    map_signif_level=T
  )  
