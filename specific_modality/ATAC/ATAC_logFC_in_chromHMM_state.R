rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(edgeR)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
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
state_num <- 15
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_summary <- data.frame()
for(tissue in tissues){
  file_dir <- paste0("result/all/ChromHMM/all_tissues_previous/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]

  if(tissue %in% c("mammarygland","ovary","uterus")){
    chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X"))),]
  }
  chromHMM_young$V2 <- chromHMM_young$V2 + 1
  chromHMM_young <- as.data.table(chromHMM_young)
  setDT(chromHMM_young)
  setkey(chromHMM_young,V1,V2,V3)
  
  search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  search_table$age <- factor(search_table$age, c("3m","24m"))
  search_table <- search_table[order(search_table$age),]
  
  tab <- read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(tab)[7:length(tab)] <-  gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
  counts <- tab[7:length(tab)]
  counts <- counts[,search_table$sample_name]
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  tab<-tab[keep,]
  tab$label <- paste0(tab$Chr,":",tab$Start,"-",tab$End)
  out <- cbind(tab[,c("label","Chr","Start","End"),drop=F],cpm(y))
  cpm_young <- out[,c("label","Chr","Start","End",search_table$sample_name[which(search_table$age=="3m")])]
  cpm_old <- out[,c("label","Chr","Start","End",search_table$sample_name[which(search_table$age=="24m")])]
  cpm_young$mean_young <- rowMeans(cpm_young[,c(5:ncol(cpm_young))])
  cpm_old$mean_old <- rowMeans(cpm_old[,c(5:ncol(cpm_old))])
  cpm_mean_summary <- merge(cpm_young[,c("label","Chr","Start","End","mean_young"),drop=F],cpm_old[,c("label","mean_old"),drop=F],by="label")
  cpm_mean_summary$logFC <- log2(cpm_mean_summary$mean_old/cpm_mean_summary$mean_young)
  cpm_mean_summary$Start <- cpm_mean_summary$Start + 1
  cpm_mean_summary <- as.data.table(cpm_mean_summary[,c("Chr","Start","End","logFC")])
  setDT(cpm_mean_summary)
  setkey(cpm_mean_summary,Chr,Start,End)
  overlaps <- as.data.frame(foverlaps(cpm_mean_summary, chromHMM_young, type = "any", nomatch = 0L))
  result <- overlaps %>%
    group_by(V4) %>%
    summarise(median_logFC = median(logFC, na.rm = TRUE))
  colnames(result)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary) == 0){
    tissue_summary <- result
  }else{
    tissue_summary <- merge(tissue_summary,result,by="V4",all=T)
  }
}
to_plot <- tissue_summary
to_plot$V4 <- factor(to_plot$V4,levels=paste0("E",1:15))
rownames(to_plot) <- to_plot$V4
to_plot <- to_plot[order(to_plot$V4),]
to_plot <- to_plot[,-1]
breaks <- c(seq(-1.5, -0.3, length.out = 40), seq(-0.29, 0.29, length.out = 20), seq(0.3,1.5, length.out = 40))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap::pheatmap(to_plot,cluster_rows = T,cluster_cols = T,show_rownames = T,breaks = breaks, color = color_palette)

for(tissue in tissues){
  ATAC_in_chromHMM_state(tissue,14)  
}
