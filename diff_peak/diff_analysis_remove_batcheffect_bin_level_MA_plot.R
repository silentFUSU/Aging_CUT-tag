rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(viridis)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K27me3"
MA_plot <- function(tissue,antibody){
  p_list <- list()
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  df <- df[,c("Geneid","logCPM","LogFC.old.young","Significant")]  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  # peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_in_young_old_merge-W1000-G3000-E100.bed"))
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_young_old_merge-W5000-G10000-E100_bedtools_filtered_500kb.bed"))
  p_list[[1]] <- ggplot(
    df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
    geom_point(aes(color = Significant),alpha=0.2, size=2) +
    scale_color_manual(values = colour) +
    labs(x="Log2(CPM)",
         y="Log2(Fold Change)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
    annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)
  df <- df[which(df$Geneid %in% peaks$V4),]
  p_list[[2]] <- ggplot(
    df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
    geom_point(aes(color = Significant),alpha=0.2,size=2) +
    scale_color_manual(values = colour) +
    labs(x="Log2(CPM)",
         y="Log2(Fold Change)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," bin in peaks"))+
    annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
    annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)
  return(p_list)
}
bin_p_list <- list()
bin_in_peaks_p_list <- list()
for(tissue in tissues){
  p_list <- MA_plot(tissue,antibody)  
  bin_p_list[[tissue]] <- p_list[[1]]
  bin_in_peaks_p_list[[tissue]] <- p_list[[2]]
}
combined_plot <- plot_a_list(bin_p_list,4,7)
ggsave(paste0("result/all/diff/",antibody,"/all_tissues_MA_plot_bin_level.png"),combined_plot,width = 35,height = 20,type="cairo")
combined_plot <- plot_a_list(bin_in_peaks_p_list,4,7)
ggsave(paste0("result/all/diff/",antibody,"/all_tissues_MA_plot_bin_young_old_merge-W5000-G10000-E100_bedtools_500kb_filtered_peaks.png"),combined_plot,width = 35,height = 20,type="cairo")

plist <- list()
for(tissue in tissues){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  df <- df[,c("Geneid","logCPM","LogFC.old.young","Significant")]
  # df <- df[which(df$Significant !="Stable"),]
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  # peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_in_young_old_merge-W1000-G3000-E100.bed"))
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_in_young_old_merge-W5000-G10000-E100_bedtools_filtered_500kb.bed"))
  df <- df[which(df$Geneid %in% peaks$V4 ),]
  # par(cex.lab = 1.5, cex.axis = 1.2)  
  plist[[tissue]] <- ggplot(df, aes(x = logCPM, y = LogFC.old.young)) +  
    geom_bin2d(bins = 200) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(  
      x = "Log2(CPM)",  
      y = "Log2(Fold Change)"
    ) + 
    scale_fill_continuous(type = "viridis", trans = "log2") + 
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," bin in peaks"))+
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")
}
combined_plot <- plot_a_list(plist,4,7)
ggsave(paste0("result/all/diff/",antibody,"/all_tissues_scatter_plot_bin_in_young_old_merge-W5000-G10000-E100_bedtools_filtered_500kb_peaks.png"),combined_plot,width = 35,height = 20,type="cairo")

