rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(edgeR)
library(data.table)

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

tissue <- "lung"
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_num <-27


condition <- "domain"
CPM_log2FC_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(condition=="domain"){
    regions <- read.table("data/samples/all/H3K27me3/bed/H3K27me3_edd_domain_merged.bed")
  }else{
    regions <- read.table(paste0("data/samples/all/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_",tissue_num,"_tissues.bed"))
  }
  tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins.counts"),header = T)
  if(tissue %in% c("mammarygland","uterus","ovary")){
    tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
    regions <- regions[which(regions$V1 %in% paste0("chr",c(1:19,"X"))),]
  }
  regions$label <- paste0(regions$V1,":",regions$V2,"-",regions$V3)
  regions <- as.data.table(regions)
  setDT(regions)
  setkey(regions,V1,V2,V3)
  
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  CPM <- as.data.frame(edgeR::cpm(counts))
  
  CPM_young <- CPM[,search_table$sample_name[which(search_table$age=="3m")]]
  CPM_old <- CPM[,search_table$sample_name[which(search_table$age=="24m")]]
  
  CPM_young$mean_young <- rowMeans(CPM_young)
  CPM_old$mean_old <- rowMeans(CPM_old)
  
  CPM_mean_summary <- merge(CPM_young[,"mean_young",drop=F],CPM_old[,"mean_old",drop=F],by="row.names")
  CPM_mean_summary$log2FC <- log2(CPM_mean_summary$mean_old/CPM_mean_summary$mean_young)
  
  CPM_mean_summary <- merge(CPM_mean_summary,tab[,c(1:4)],by.x="Row.names",by.y="Geneid")
  CPM_mean_summary <- as.data.table(CPM_mean_summary)
  setDT(CPM_mean_summary)
  setkey(CPM_mean_summary,Chr,Start,End)
  
  overlaps <- as.data.frame(foverlaps(regions, CPM_mean_summary, type = "any", nomatch = 0L))
  # CPM_mean_summary <- overlaps %>%
  #   group_by(label) %>%
  #   summarise(
  #     median_young = median(mean_young, na.rm = TRUE),
  #     median_old = median(mean_old, na.rm = TRUE)
  #   )
  CPM_mean_summary <- overlaps %>%
    group_by(label) %>%
    summarise(
      log2FC = median(log2FC, na.rm = TRUE)
    )
  colnames(CPM_mean_summary)[which(colnames(CPM_mean_summary)=="log2FC")] <- tissue_label_change(tissue)
  colnames(CPM_mean_summary)[1] <- "Geneid"
  if(nrow(CPM_log2FC_summary )==0){
    CPM_log2FC_summary <- CPM_mean_summary[,c("Geneid",tissue_label_change(tissue))]
  }else{
    CPM_log2FC_summary <- merge(CPM_log2FC_summary,CPM_mean_summary[,c("Geneid",tissue_label_change(tissue))],by="Geneid",all=T)
  }
}
rownames(CPM_log2FC_summary) <- CPM_log2FC_summary$Geneid
CPM_log2FC_summary <- CPM_log2FC_summary[,-1]
CPM_log2FC_summary$row_mean <- rowMeans(CPM_log2FC_summary, na.rm = TRUE)
col_means <- colMeans(as.matrix(CPM_log2FC_summary[, -ncol(CPM_log2FC_summary)]), na.rm = TRUE)
CPM_log2FC_summary_sorted_rows <- CPM_log2FC_summary %>%
  arrange(row_mean)
CPM_log2FC_summary_sorted_rows <- CPM_log2FC_summary_sorted_rows[ , -ncol(CPM_log2FC_summary_sorted_rows)]
CPM_log2FC_summary_sorted <- CPM_log2FC_summary_sorted_rows[, order(col_means)]
CPM_log2FC_summary <- CPM_log2FC_summary_sorted

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(-5, -0.51, length.out = 40), seq(-0.5, 0.5, length.out = 20), seq(0.51, 5, length.out = 40))
p <- pheatmap::pheatmap(CPM_log2FC_summary,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")
p <- pheatmap::pheatmap(CPM_log2FC_summary,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks = breaks, color = color_palette, clustering_distance_cols="manhattan",clustering_distance_rows="manhattan")
