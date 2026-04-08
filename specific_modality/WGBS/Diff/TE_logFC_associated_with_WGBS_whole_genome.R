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

tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 

tissue_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_200kb_bins_cross_comparison.csv",row.names = 1)
H3K9me3_tissue_order_label <- c() 
annotation_col <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- t_search_table$sample_name[which(t_search_table$age=="3M")]
  old_samples <- t_search_table$sample_name[which(t_search_table$age=="24M")]
  combinations <- as.data.frame(expand.grid(young = young_samples, old = old_samples))
  if(tissue=="bonemarrow"){
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Bone.Marrow",".",paste0(combinations$old,".",combinations$young)))  
  }else if(tissue=="mammarygland"){
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0("Mammary.Gland",".",paste0(combinations$old,".",combinations$young)))  
  }
  else{
    H3K9me3_tissue_order_label <- c(H3K9me3_tissue_order_label, paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
    t_annotation_col <- data.frame(tissue=tissue_label_change(tissue),sample=paste0(tissue_label_change(tissue),".",paste0(combinations$old,".",combinations$young)))  
  }
  annotation_col <- rbind(annotation_col,t_annotation_col)
}
rownames(tissue_summary) <- tissue_summary$label

tissue_mean_summary <- data.frame() 
annotation <- annotation_col
for(tissue in tissues){
  t_annotation <- annotation[which(annotation$tissue==tissue_label_change(tissue)),]
  t_tissue_mean_summary <- tissue_summary[,t_annotation$sample]
  t_tissue_mean_summary$mean_delta <- rowMeans(t_tissue_mean_summary) 
  t_tissue_mean_summary$label <- rownames(t_tissue_mean_summary)
  t_tissue_mean_summary <- t_tissue_mean_summary[,c("label","mean_delta")]
  colnames(t_tissue_mean_summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_mean_summary)==0){
    tissue_mean_summary <- t_tissue_mean_summary  
  }else{
    tissue_mean_summary <- merge(tissue_mean_summary,t_tissue_mean_summary,by="label")
  }
}

to_plot_WGBS_order_long <- reshape2::melt(tissue_mean_summary)

medians <- to_plot_WGBS_order_long %>%
  group_by(variable) %>%
  summarise(median_value = median(value, na.rm = TRUE))
medians <- medians[-which(medians$variable %in% c("Mammary Gland","Uterus","Ovary")),]
medians <- medians[order(medians$median_value),]
medians$rank <- c(1:nrow(medians))
medians <- as.data.frame(medians)

gene <- "ERVK"
summary <- data.frame()
for(tissue in tissues){
  df <- read.delim(paste0("data/samples/RNA/",tissue,"/TEcount/combined.cntTable"),row.names = 1)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SRR[0-9]+|HM[0-9]+).*"
  colnames(df) <- gsub(pattern, "\\1", colnames(df))
  CPM <- as.data.frame(cpm(df))
  matching_rows <- grep(paste0(":",gene,":"), rownames(CPM))
  CPM <- CPM[matching_rows, ]
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$sample_name %in% colnames(df)),]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(df))
  search_table <- search_table[order(search_table$sample_name),]
  young_samples <- CPM[,search_table$sample_name[which(search_table$age=="3m")]]
  young_samples <- colSums(young_samples)
  young_samples <- mean(young_samples)
  old_samples <- CPM[,search_table$sample_name[which(search_table$age=="24m")]]
  old_samples <- colSums(old_samples)
  old_samples <- mean(old_samples)
  t_summary <- data.frame(tissue=tissue_label_change(tissue),logFC=log2(old_samples/young_samples))
  summary <- rbind(summary,t_summary)
}
summary <- summary[order(summary$logFC),]
summary$tissue <- factor(summary$tissue,levels = summary$tissue)

colnames(summary)[2] <- "logFC_TE" 
colnames(medians)[1:2] <-c("tissue","delta_WGBS" ) 
to_plot <- merge(summary,medians,by="tissue")

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(to_plot$tissue))
cor.test(to_plot$logFC_TE,to_plot$delta_WGBS,method="spearman")
ggplot(to_plot,aes(x=delta_WGBS,y=logFC_TE,color =tissue))+    
  geom_jitter(size = 3, alpha = 0.7)+
  scale_color_manual(values = color) +
  ggtitle(paste0(gene," expression change relationship with DNA methylation change"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("DNA methylation Delta")+labs(fill = "", color = "") +ylab(paste0(gene," log2(Fold change)"))+
  scale_x_reverse()

