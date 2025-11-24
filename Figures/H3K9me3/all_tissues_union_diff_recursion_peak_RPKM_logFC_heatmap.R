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

# antibody <- "H3K9me3"
common_increase <- read.csv(paste0("data/samples/all/H3K9me3/common_increase-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_increase <- common_increase[which((common_increase$end - common_increase$start +1) > 200000),]
common_decrease <- read.csv(paste0("data/samples/all/H3K9me3/common_decrease-W5000-G10000-E100_recursion_union_peaks_after_remove_batch_effect.csv"))
common_decrease <- common_decrease[which((common_decrease$end - common_decrease$start +1) > 200000),]
regions <- unique(c(common_increase$Geneid[which(common_increase$n>0)],common_decrease$Geneid[which(common_decrease$n>0)]))
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "ATAC"
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
  colnames(rpkm_mean_summary)[which(colnames(rpkm_mean_summary)=="log2FC")] <- tissue_label_change(tissue)
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  if(nrow(diff_summary )==0){
    diff_summary <- rpkm_mean_summary[,c(1,4)]
  }else{
    diff_summary <- merge(diff_summary,rpkm_mean_summary[,c(1,4)],by="Geneid",all=T)
  }
}

rownames(diff_summary) <- diff_summary$Geneid
diff_summary <- diff_summary[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
# p <- pheatmap::pheatmap(diff_summary,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")

kmeans_result <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
diff_summary <- merge(diff_summary,kmeans_result,by.x="row.names",by.y="X")
rownames(diff_summary) <- diff_summary$Row.names
diff_summary <- diff_summary[order(diff_summary$cluster),]
to_plot <- diff_summary[,-c(1,ncol(diff_summary))]

annotation <- kmeans_result
rownames(annotation) <- annotation$X
annotation <- annotation[,-1,drop=F]
annotation$cluster <- as.character(annotation$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissue_order]
pheatmap::pheatmap(to_plot,cluster_rows = F,main = antibody,cluster_cols = F,show_rownames = F,breaks = breaks, annotation_row = annotation,color = color_palette, clustering_distance_cols="manhattan")
