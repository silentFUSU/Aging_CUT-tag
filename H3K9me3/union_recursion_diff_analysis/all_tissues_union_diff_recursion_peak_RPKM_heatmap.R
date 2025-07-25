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

variance <- list()
antibody <- "H3K9me3"
common_increase <- read.csv(paste0("data/samples/all/",antibody,"/common_increase-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_increase <- common_increase[which((common_increase$end - common_increase$start +1) > 200000),]
common_decrease <- read.csv(paste0("data/samples/all/",antibody,"/common_decrease-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_decrease <- common_decrease[which((common_decrease$end - common_decrease$start +1) > 200000),]
regions <- unique(c(common_increase$Geneid[which(common_increase$n>0)],common_decrease$Geneid[which(common_decrease$n>0)]))
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

window_size="5000"
gap_size="10000"
age <- "3m"
rpkm_summary <- data.frame()

for(tissue in tissues){
  df <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion.counts"),skip=1)
  df <- df[which(df$Geneid %in% regions),]
  if(tissue %in% c("mammarygland","ovary","uterus")){
    df <- df[-which(df$Chr == "chrY"),]
  }
  
  counts = df[,c(7:ncol(df))]
  rownames(counts)= df$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  
  counts <- counts[,search_table$sample_name]
  total_reads <- colSums(counts)
  feature_lengths <- df$Length
  rpkm <- counts
  for (i in seq_along(total_reads)) {
    rpkm[, i] <- rpkm[, i] / (df$Length / 1000) / (total_reads[i] / 1e6)
  }
  
  rownames(rpkm) <- rownames(counts)
  colnames(rpkm) <- colnames(counts)
  
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody=="H3K9me3" & search_table$age==age),]
  rpkm <- rpkm[,c(search_table$sample_name)]
  
  rpkm$rpkm <- rowSums(rpkm)/nrow(search_table)
  rpkm$Geneid <- rownames(rpkm)
  rpkm <- rpkm[,c("Geneid","rpkm")]
  colnames(rpkm)[2] <- tissue_label_change(tissue)
  
  if(nrow(rpkm_summary) == 0){
    rpkm_summary <- rpkm
  }else{
    rpkm_summary <- merge(rpkm_summary,rpkm,by="Geneid",all=T)    
  }
}
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
annotation <- annotation[order(annotation$cluster),]
rownames(rpkm_summary) <- rpkm_summary$Geneid
rpkm_summary <- rpkm_summary[,-1]
rpkm_summary <- rpkm_summary[annotation$X,]
tissues_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung","Mammary Gland",
                   "Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
rpkm_summary <- rpkm_summary[annotation$X,tissues_order]
rownames(annotation) <- annotation$X
annotation <- annotation[,-1,drop=F]
rpkm_summary <- log2(rpkm_summary)
annotation$cluster <- as.character(annotation$cluster)
breaks <- c(seq(-6, -2.1, length.out = 40), seq(-2, 2, length.out = 20), seq(2.1, 6, length.out = 40))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
pheatmap::pheatmap(rpkm_summary,cluster_rows = F,cluster_cols = F,main=age,annotation_row = annotation,breaks = breaks,show_rownames = F,color = color_palette)
breaks <- c(seq(-4, -2.1, length.out = 40), seq(-2, 2, length.out = 20), seq(2.1, 4, length.out = 40))
pheatmap::pheatmap(rpkm_summary,cluster_rows = F,cluster_cols = F,main=age,breaks=breaks,annotation_row = annotation,scale="row",show_rownames = F,color = color_palette)

variance[[age]] <-  apply(rpkm_summary, 1, var)
variance_summary <- data.frame(Geneid=rownames(annotation),Young=variance$`3m`,Old=variance$`24m`)
variance_summary <- variance_summary[,-1]
breaks <- c(seq(-1, -0.6, length.out = 40), seq(-0.5, 0.5, length.out = 20), seq(0.6, 1, length.out = 40))
pheatmap::pheatmap(variance_summary,cluster_rows = F,cluster_cols = F,main="Variation between tissues ",breaks=breaks,annotation_row = annotation,scale="row",show_rownames = F,color = color_palette)
to_plot <- variance_summary
to_plot$Geneid <- rownames(to_plot)
to_plot <- reshape2::melt(to_plot)
to_plot$variable <- factor(to_plot$variable,levels=c("Young","Old")) 
ggplot(to_plot, aes(x = variable, y = value,fill=variable))+
  geom_boxplot(outliers = F) +
  labs(x = "Age", y = "Variation") +
  theme_minimal()+
  theme(
    text = element_text(size = 14),        
    axis.title = element_text(size = 16),  
    axis.text = element_text(size = 12),   
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)  
  )+
  ggtitle("Variation between tissues")+
  geom_signif(
    comparisons = list(c("Young", "Old")),
    map_signif_level = F,
    textsize = 4,test = "wilcox.test",
    y_position = c(0.3),
    tip_length = c(1/100)
  )
