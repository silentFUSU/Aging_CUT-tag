rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(edgeR)
library(data.table)
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

H3K9me3_peaks_CPM_relationship_H3K27me3 <- function(tissue){
  tab<- read.delim(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_merge-W1000-G3000-E100_compress.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$age=="3m" & search_table$antibody=="H3K9me3"),]
  
  counts <- counts[,which(colnames(counts) %in% search_table$sample_name)]
  cpm <- as.data.frame(cpm(counts,log = T))
  cpm$mean <- (cpm[,1]+cpm[,2])/2
  H3K9me3 <- cbind(tab[,c(1:6)],enrichment=cpm$mean)
  # peak_lengths <- tab$Length
  # peak_lengths <- peak_lengths/1000
  # total_reads <- colSums(counts)
  # total_reads <- total_reads/ 1e6
  # rpkm_matrix <- apply(counts, 2, function(x) {  
  #   x / (peak_lengths * total_reads)  
  # }) 
  # rpkm <- as.data.frame(rpkm_matrix)
  # rpkm$mean <- (rpkm[,1]+rpkm[,2])/2
  # H3K9me3 <- cbind(tab[,c(1:6)],enrichment=rpkm$mean)
  
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K27me3 <- H3K27me3[,c("Chr","Start","End","LogFC.old.young","Significant")]
  
  H3K9me3$Start <- H3K9me3$Start + 1
  H3K9me3 <- as.data.table(H3K9me3)
  setDT(H3K9me3)
  setkey(H3K9me3,Chr,Start,End)
  
  H3K27me3$Start <- H3K27me3$Start + 1
  H3K27me3 <- as.data.table(H3K27me3)
  setDT(H3K27me3)
  setkey(H3K27me3,Chr,Start,End)
  
  overlaps <- foverlaps(H3K9me3,H3K27me3, type = "any", nomatch = 0L)
  overlaps <- as.data.frame(overlaps)
  to_plot<- overlaps %>%  
    group_by(Geneid) %>%  
    summarize(mean_LogFC = mean(LogFC.old.young, na.rm = TRUE))  
  to_plot <- merge(to_plot,H3K9me3[,c("Geneid","enrichment")],by="Geneid")
  
  
  result <- overlaps %>%  
    count(Geneid, Significant) %>%  
    group_by(Geneid) %>%  
    mutate(proportion = n / sum(n)) %>%  
    arrange(Geneid, desc(proportion), factor(Significant, levels = c("Up", "Down", "Stable"))) %>%  
    slice(1) %>%  
    select(Geneid, Significant)  
  to_plot <- merge(to_plot,result,by="Geneid")
  to_plot <- merge(to_plot,tab[,c("Geneid","Chr","Start","End","Length")], by="Geneid")
  up_count <- sum(to_plot$Significant == "Up")  
  down_count <- sum(to_plot$Significant == "Down")  
  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(to_plot , aes(x = enrichment, y =mean_LogFC, color=Significant)) +  
    geom_point() +  
    scale_color_manual(values = colour)+
    theme_bw()+
    ggtitle(tissue_label_change(tissue),"H3K27me3 log2(FC) ~ H3K9me3 peaks enrichment") +
    labs(x = "H3K9me3 enrichment", y = "log2(H3K27me3 Fold Change)") +
    theme(  
      plot.title = element_text(size = 14),     
      axis.title.x = element_text(size = 12),   
      axis.title.y = element_text(size = 12),   
      axis.text.x = element_text(size = 10),   
      axis.text.y = element_text(size = 10)    
    ) +
    annotate("text", x = Inf, y = Inf,  
             label = paste("Up:", up_count),  
             hjust = 1.5, vjust = 1.5, size = 5, color = "black") +  
    annotate("text", x = Inf, y = -Inf,  
             label = paste("Down:", down_count),  
             hjust = 1.5, vjust = -1.5, size = 5, color = "black")
  return(p)
}
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- H3K9me3_peaks_CPM_relationship_H3K27me3(tissue)  
}
combined_plot <- plot_a_list(p_list,no_of_cols = 7,no_of_rows = ceiling(length(p_list)/7))
ggsave(paste0("result/all/H3K27me3_H3K9me3/all_tissues_H3K27me3_log2FC_H3K9me3_peaks_CPM.png"),width = 35,height = ceiling(length(p_list)/7)*5, type="cairo")

H3K9me3_peaks_length_relationship_H3K27me3 <- function(tissue){
  tab<- read.delim(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_merge-W1000-G3000-E100.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$age=="3m" & search_table$antibody=="H3K9me3"),]
  
  H3K9me3 <- cbind(tab[,c(1:6)])
  # peak_lengths <- tab$Length
  # peak_lengths <- peak_lengths/1000
  # total_reads <- colSums(counts)
  # total_reads <- total_reads/ 1e6
  # rpkm_matrix <- apply(counts, 2, function(x) {  
  #   x / (peak_lengths * total_reads)  
  # }) 
  # rpkm <- as.data.frame(rpkm_matrix)
  # rpkm$mean <- (rpkm[,1]+rpkm[,2])/2
  # H3K9me3 <- cbind(tab[,c(1:6)],enrichment=rpkm$mean)
  
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K27me3 <- H3K27me3[,c("Chr","Start","End","LogFC.old.young","Significant")]
  
  H3K9me3$Start <- H3K9me3$Start + 1
  H3K9me3 <- as.data.table(H3K9me3)
  setDT(H3K9me3)
  setkey(H3K9me3,Chr,Start,End)
  
  H3K27me3$Start <- H3K27me3$Start + 1
  H3K27me3 <- as.data.table(H3K27me3)
  setDT(H3K27me3)
  setkey(H3K27me3,Chr,Start,End)
  
  overlaps <- foverlaps(H3K9me3,H3K27me3, type = "any", nomatch = 0L)
  overlaps <- as.data.frame(overlaps)
  to_plot<- overlaps %>%  
    group_by(Geneid) %>%  
    summarize(mean_LogFC = mean(LogFC.old.young, na.rm = TRUE))  
  to_plot <- merge(to_plot,H3K9me3[,c("Geneid","Length")],by="Geneid")
  
  
  result <- overlaps %>%  
    count(Geneid, Significant) %>%  
    group_by(Geneid) %>%  
    mutate(proportion = n / sum(n)) %>%  
    arrange(Geneid, desc(proportion), factor(Significant, levels = c("Up", "Down", "Stable"))) %>%  
    slice(1) %>%  
    select(Geneid, Significant)  
  to_plot <- merge(to_plot,result,by="Geneid")
  up_count <- sum(to_plot$Significant == "Up")  
  down_count <- sum(to_plot$Significant == "Down")  
  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(to_plot , aes(x = log2(Length), y =mean_LogFC, color=Significant)) +  
    geom_point() +  
    scale_color_manual(values = colour)+
    theme_bw()+
    ggtitle(tissue_label_change(tissue),"H3K27me3 log2(FC) ~ H3K9me3 peaks log2(length)") +
    labs(x = "H3K9me3 Length", y = "log2(H3K27me3 Fold Change)") +
    theme(  
      plot.title = element_text(size = 14),     
      axis.title.x = element_text(size = 12),   
      axis.title.y = element_text(size = 12),   
      axis.text.x = element_text(size = 10),   
      axis.text.y = element_text(size = 10)    
    ) +
    annotate("text", x = Inf, y = Inf,  
             label = paste("Up:", up_count),  
             hjust = 1.5, vjust = 1.5, size = 5, color = "black") +  
    annotate("text", x = Inf, y = -Inf,  
             label = paste("Down:", down_count),  
             hjust = 1.5, vjust = -1.5, size = 5, color = "black")
  return(p)
}
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- H3K9me3_peaks_length_relationship_H3K27me3(tissue)  
}
combined_plot <- plot_a_list(p_list,no_of_cols = 7,no_of_rows = ceiling(length(p_list)/7))
ggsave(paste0("result/all/H3K27me3_H3K9me3/all_tissues_H3K27me3_log2FC_H3K9me3_peaks_length.png"),width = 35,height = ceiling(length(p_list)/7)*5, type="cairo")


H3K9me3_bin_in_peaks_CPM_relationship_H3K27me3 <- function(tissue){
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  colnames(H3K9me3)[which(colnames(H3K9me3)=="logCPM")] <- "log2(H3K9me3_CPM)"
  
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  colnames(H3K27me3)[which(colnames(H3K27me3)=="LogFC.old.young")] <- "log2(H3K27me3_Fold_Change)"
  
  H3K9me3_peak <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  H3K9me3 <- H3K9me3[which(H3K9me3$Geneid %in% H3K9me3_peak$V4),]
 
  df <- merge(H3K9me3[,c("Geneid","log2(H3K9me3_CPM)")],H3K27me3[,c("Geneid","log2(H3K27me3_Fold_Change)","Significant")], by="Geneid")
  up_count <- sum(df$Significant == "Up")  
  down_count <- sum(df$Significant == "Down")  
  color <- setNames(c("red","grey","blue"),c("Up","Stable","Down"))
  p <- ggplot(df, aes(x = `log2(H3K9me3_CPM)`, y = `log2(H3K27me3_Fold_Change)`,color=Significant)) +  
    geom_point() +  
    scale_color_manual(values = color)+
    theme_bw()+
    ggtitle(tissue_label_change(tissue),"H3K27me3 log2(FC) ~ H3K9me3 CPM in H3K9me3 peaks") +
    labs(x = "log2(H3K9me3 CPM)", y = "log2(H3K27me3 Fold Change)") +
    theme(  
      plot.title = element_text(size = 14),     
      axis.title.x = element_text(size = 12),   
      axis.title.y = element_text(size = 12),   
      axis.text.x = element_text(size = 10),   
      axis.text.y = element_text(size = 10)    
    ) +
    annotate("text", x = Inf, y = Inf,  
             label = paste("Up:", up_count),  
             hjust = 1.5, vjust = 1.5, size = 5, color = "black") +  
    annotate("text", x = Inf, y = -Inf,  
             label = paste("Down:", down_count),  
             hjust = 1.5, vjust = -1.5, size = 5, color = "black")
  return(p)
}
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- H3K9me3_bin_in_peaks_CPM_relationship_H3K27me3 (tissue)
}

combined_plot <- plot_a_list(p_list,no_of_cols = 7,no_of_rows = ceiling(length(p_list)/7))
ggsave(paste0("result/all/H3K27me3_H3K9me3/all_tissues_H3K27me3_log2FC_H3K9me3_CPM_in_peaks.png"),width = 35,height = ceiling(length(p_list)/7)*5, type="cairo")

