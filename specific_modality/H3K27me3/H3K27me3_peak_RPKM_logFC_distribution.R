rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(enrichplot)
options(scipen = 0) 
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
antibody <- "H3K27me3"
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver","lung","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT","muscle")
rpkm_log2FC_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_old_merge-W5000-G10000-E100.counts"),header = T)
  tab <- tab[which(tab$Length < 50000),]
  summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_old_merge-W5000-G10000-E100.counts.summary"),header = T,row.names = 1)
  if(tissue %in% c("mammarygland","uterus","ovary")){
    tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
  }
  
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  colnames(summary) <- gsub(pattern,"\\1",colnames(summary))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  summary <- summary[-2,search_table$sample_name]
  total_reads <- colSums(summary)
  length <- as.numeric(tab$Length)
  
  rpkm <- sweep(counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  
  rpkm_young <- rpkm[,search_table$sample_name[which(search_table$age=="3m")]]
  rpkm_old <- rpkm[,search_table$sample_name[which(search_table$age=="24m")]]
  
  rpkm_young$mean_young <- rowMeans(rpkm_young)
  rpkm_old$mean_old <- rowMeans(rpkm_old)
  
  rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
  rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
  
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_old_merge-W5000-G10000-E100_diff_after_remove_batch_effect.csv"))
  rpkm_mean_summary <- rpkm_mean_summary[which(rpkm_mean_summary$Row.names %in% df$Geneid[which(df$Significant=="Down")]),]
  rpkm_mean_summary$tissue <- tissue_label_change(tissue)
  t_rpkm_log2FC_summary <- rpkm_mean_summary[,c("tissue","log2FC")]
  rpkm_log2FC_summary <- rbind(rpkm_log2FC_summary,t_rpkm_log2FC_summary)
}
to_plot <- rpkm_log2FC_summary
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))
average_log2FC <- to_plot %>%
  group_by(tissue) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE))
average_log2FC <- average_log2FC[order(average_log2FC$mean_log2FC),]
to_plot$tissue <- factor(to_plot$tissue,levels=average_log2FC$tissue)
ggplot(to_plot, aes(x = tissue, y = log2FC,fill=tissue)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = color)+
  labs(title = "Boxplot of Length by Tissue",
       x = "Tissue",
       y = "log2FC") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))


