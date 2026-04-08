rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(factoextra)
library(cluster)
library(umap)
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

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K27me3"
diff_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(antibody == "H3K9me3"){
    tab <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.counts.summary"),header = T,row.names = 1)
  }else if(antibody == "RNA"){
    search_table <- read.csv("data/samples/all/RNA_search_table.csv")
    tab <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_H3K9me3_peaks.counts"),header = T) 
    summary <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)  
  }else if(antibody == "ATAC"){
    search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
    tab <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/",tissue,"_H3K9me3_peaks.counts"),header = T) 
    summary <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/",tissue,"_H3K9me3_peaks.counts.summary"),header=T)
  }else{
    tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts"),header = T)
    summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts.summary"),header = T)
  }
  if(tissue %in% c("mammarygland","uterus","ovary")){
    tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
  }
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
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
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  H3K9me3_diff <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  H3K9me3_diff <- H3K9me3_diff[,c("Geneid","Length","LogFC.old.young","Significant")]
  colnames(H3K9me3_diff)[4] <- "histone_significant"
  rpkm_mean_summary <- merge(H3K9me3_diff,rpkm_mean_summary,by="Geneid")
  colnames(rpkm_mean_summary)[c(3,7)] <- c("LogFC.old.young","logFC")
  rpkm_mean_summary$tissue <- tissue_label_change(tissue)
  diff_summary <- rbind(diff_summary,rpkm_mean_summary)
}

diff_summary <- diff_summary[which(diff_summary$Length > 200000),]
to_plot <- diff_summary %>%  
  group_by(tissue, histone_significant) %>%  
  summarise(  
    median_logFC = median(logFC, na.rm = TRUE),  
    count = n()  
  ) %>%   
  mutate(median_logFC = ifelse(count < 10, NA, median_logFC))

to_plot <- to_plot[,-ncol(to_plot)]
to_plot <- as.data.frame(to_plot)
to_plot <- reshape2::dcast(to_plot, tissue ~ histone_significant, value.var = "median_logFC")  
rownames(to_plot) <- to_plot$tissue
to_plot <- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-0.5, -0.21, length.out = 40), seq(-0.2, 0.2, length.out = 20), seq(0.21, 0.5, length.out = 40))  
tissues_order <- c("Testis","Tongue","Stomach","Cecum","Colon","Pancreas","Ileum","Liver","Heart","Jejunum","Hippocampus","Muscle","Bone Marrow","Ovary","iWAT","Cortex","Aorta","Uterus","Bladder","Spleen","Thymus","Kidney","Skin","Cerebellum","Lung","BAT","Mammary Gland")
to_plot <- to_plot[tissues_order,]
to_plot <- to_plot[,c("Up","Stable","Down")]
pheatmap::pheatmap(to_plot,breaks = breaks,color = color_palette,cluster_rows = F,cluster_cols = F,na_col = "grey")

p_value_summary <- data.frame(
  Up = rep(NA, 27), 
  Stable = rep(NA, 27), 
  Down = rep(NA, 27)  
)
rownames(p_value_summary) <- sapply(tissues, tissue_label_change)
for(tissue in tissues){
  t_Up_summary <-  diff_summary[which(diff_summary$tissue == tissue_label_change(tissue) & diff_summary$histone_significant=="Up"),]
  t_Down_summary <- diff_summary[which(diff_summary$tissue == tissue_label_change(tissue) & diff_summary$histone_significant=="Down"),]
  t_Stable_summary <- diff_summary[which(diff_summary$tissue == tissue_label_change(tissue) & diff_summary$histone_significant=="Stable"),]
  if(nrow(t_Up_summary) > 10){
    test <- wilcox.test(t_Up_summary$logFC,t_Stable_summary$logFC)
    p_value_summary[tissue_label_change(tissue),"Up"] <- test$p.value
  }
  if(nrow(t_Down_summary) > 10){
    test <- wilcox.test(t_Down_summary$logFC,t_Stable_summary$logFC)
    p_value_summary[tissue_label_change(tissue),"Down"] <- test$p.value
  }
}

mark_significance <- function(p_value) {
  if (is.na(p_value)) {
    return(NA)
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return(NA)
  }
}
p_value_summary <- p_value_summary %>%
  mutate(
    Up = sapply(Up, mark_significance),
    Stable = sapply(Stable, mark_significance),
    Down = sapply(Down, mark_significance)
  )
p_value_summary <- p_value_summary[tissues_order,]
to_plot$tissue <- rownames(to_plot)
df_long <- reshape2::melt(to_plot)
names(df_long) <- c("Tissue", "Type", "Value")

p_value_summary$tissue <- rownames(p_value_summary)
p_value_long <- reshape2::melt(p_value_summary,id.vars = "tissue")
names(p_value_long) <- c("Tissue", "Type", "Label")
merged_data <- merge(df_long, p_value_long, by = c("Tissue", "Type"), all.x = TRUE)
merged_data$Tissue <- factor(merged_data$Tissue,levels=rev(tissues_order))
merged_data$Value[which(merged_data$Value > 1)] <- 1
merged_data$Value[which(merged_data$Value < -1)] <- -1

ggplot(merged_data, aes(x = Type, y = Tissue, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-1, 1), midpoint = 0) +
  theme_minimal() +
  xlab("H3K9me3 change condition")+
  ggtitle(paste0(antibody," log2(Fold change)"))+
  geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

