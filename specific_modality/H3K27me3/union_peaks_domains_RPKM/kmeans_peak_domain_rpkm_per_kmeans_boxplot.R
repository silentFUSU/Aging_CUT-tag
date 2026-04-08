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

condition <- "domain"
tissue_num <- 27
if(condition=="domain"){
  regions <- read.csv("data/samples/all/H3K27me3/edd_domain_merged/kmeans_annotation_RPKM.csv")
}else{
  regions <- read.csv(paste0("data/samples/all/H3K27me3/peaks_merged/kmeans_annotation_larger_",tissue_num,"_tissues_RPKM.csv"))
}
colnames(regions)[1] <- "label"
antibody <- "H3K27me3"
window_size=5000
gap_size=10000
summary <- data.frame()
age="24m"
for(tissue in tissues){
  if(condition == "domain"){
    df <- read.delim(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts"),skip=1)
    df_summary <- read.delim(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts.summary"),header = T)
  }else{
    df <- read.delim(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.counts"),skip=1)
    df_summary <- read.delim(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.counts.summary"),header = T)
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
  df <- merge(df,regions[,c("label","cluster")],by="label")
  df$tissue <- tissue_label_change(tissue)
  # df$cluster[which(df$Chr=="chrY" & df$cluster=="kmeans1")] <- "chrY_kmeans1"
  summary <- rbind(summary,df)
}

tissue_order <- c("Lung","Cerebellum","BAT","Muscle","Heart","Aorta","Skin","Kidney","Hippocampus","Cortex","Liver","Tongue","Uterus","Testis","Bladder","Ovary","Colon","Stomach",
                  "Thymus","Cecum","Jejunum","Pancreas","Bone Marrow","Ileum","Spleen","iWAT","Mammary Gland")
summary$tissue <- factor(summary$tissue,levels=tissue_order)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
summary$cluster <- paste0("kmeans",summary$cluster)
summary$cluster <- factor(summary$cluster,levels=c("kmeans1","kmeans2","kmeans3"))
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
  )+ylim(0,1.5)

