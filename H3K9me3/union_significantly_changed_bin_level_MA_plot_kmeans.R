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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
antibody <- "H3K9me3"
kmeans <- "1"
MA_plot <- function(tissue,antibody,kmeans){
  regions <- read.csv("data/samples/all/H3K9me3/kmeans_H3K9me3_union_significantly_changed_bins.csv")
  regions <- regions[which(regions$cluster==kmeans),]
  p_list <- list()
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  df <- df[,c("Geneid","logCPM","LogFC.old.young","Significant")]  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  ymax <- max(abs(df$LogFC.old.young))
  ymax <- ymax + 0.5
  df <- df[which(df$Geneid %in% regions$X),]
  p <- ggplot(
    df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
    geom_point(aes(color = Significant),alpha=0.2,size=2) +
    scale_color_manual(values = colour) +
    labs(x="Log2(CPM)",
         y="Log2(Fold Change)") +
    theme_bw()+
    ylim(-ymax,ymax) +
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," bin kmeans",kmeans))+
    annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
    annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)+
    geom_hline(yintercept = 0, linetype = "dashed")  
  return(p)
}

tissues<- c("mammarygland","ovary","ileum","kidney","spleen","CB","cecum","colon","jejunum","thymus","bonemarrow",
            "pancreas","iWAT","liver","bladder","testis","aorta","brain","tongue","Hip","muscle","uterus","heart",
            "stomach","skin","BAT","lung")
for(kmeans in c(1:6)){
  p_list <- list()
  for(tissue in tissues){
    p_list[[tissue]] <- MA_plot(tissue,antibody,kmeans)  
  }
  combined_plot <- plot_a_list(p_list,4,7)
  ggsave(paste0("result/all/diff/",antibody,"/all_tissues_MA_plot_bin_level_union_changed_kmeans",kmeans,".png"),combined_plot,width = 35,height = 20,type="cairo")
}
