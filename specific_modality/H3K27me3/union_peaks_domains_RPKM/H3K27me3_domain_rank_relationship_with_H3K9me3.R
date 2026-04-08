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
H3K9me3_rank <- data.frame(
  tissue = c("Lung", "Cerebellum", "BAT", "Muscle", "Heart", "Aorta", "Skin", 
             "Kidney", "Hippocampus", "Cortex", "Liver", "Tongue", "Testis", 
             "Bladder", "Pancreas", "Cecum", "Spleen", "Stomach", "Colon", 
             "Bone Marrow", "Jejunum", "iWAT", "Thymus", "Ileum"),
  Median_Log2_FC = c(-0.36, -0.33, -0.32, -0.19, -0.19, -0.18, -0.18, -0.17, -0.15,
                     -0.10, -0.10, -0.09, -0.08, -0.06, -0.03, -0.03, -0.02, 0.01, 
                     0.02, 0.03, 0.04, 0.05, 0.05, 0.09),
  H3K9me3_rank = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
           21, 22, 23, 24)
)

antibody <- "H3K27me3"
annotation <- read.csv("data/samples/all/H3K27me3/edd_domain_merged/kmeans_annotation_RPKM.csv")
regions <- annotation[which(annotation$cluster=="2"),]
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver","lung","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT","muscle")
diff_summary <- data.frame()
rpkm_log2FC_summary <- data.frame()
condition <- "domain"
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(condition=="domain"){
    tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_edd_domain_merged.counts.summary"),header = T,row.names = 1)
  }else{
    tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_20_tissues.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_20_tissues.counts.summary"),header = T,row.names = 1)
  }
  
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
  colnames(rpkm_mean_summary)[which(colnames(rpkm_mean_summary)=="log2FC")] <- tissue_label_change(tissue)
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  rpkm_mean_summary <- rpkm_mean_summary[which(rpkm_mean_summary$Geneid %in% regions$X),]
  if(nrow(rpkm_log2FC_summary )==0){
    rpkm_log2FC_summary <- rpkm_mean_summary[,c(1,4)]
  }else{
    rpkm_log2FC_summary <- merge(rpkm_log2FC_summary,rpkm_mean_summary[,c(1,4)],by="Geneid",all=T)
  }
}
diff_summary <- rpkm_log2FC_summary
rownames(diff_summary) <- diff_summary$Geneid
diff_summary <- diff_summary[,-1]

medians <- apply(diff_summary, 2, median, na.rm = TRUE)
medians <- data.frame(tissues=colnames(diff_summary),median=medians)
medians <- medians[order(medians$median,decreasing = T),]
medians$rank <- 1:nrow(medians)
colnames(medians) <- c("tissue","logFC","H3K27me3_rank")

to_plot <- merge(H3K9me3_rank,medians,by="tissue")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))
cor.test(to_plot$H3K27me3_rank,to_plot$H3K9me3_rank)
ggplot(to_plot,aes(x=H3K27me3_rank,y=H3K9me3_rank,color = tissue))+  
  geom_text(aes(label = tissue),color = "black", vjust = 0, hjust = 0, size = 4) +
  geom_jitter( size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle("H3K27me3 domain kmeans2 rank relationship with H3K9me3 peaks kmeans1 rank")+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("H3K27me3 rank")+ylab("H3K9me3 rank")+labs(fill = "", color = "")+
  scale_x_reverse() +
  scale_y_reverse()

