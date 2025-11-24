rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
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
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }
  }
  return(tissue_label)
}
bin_size <- function(antibody){
  if(antibody %in% c("H3K27ac","H3K4me3","H3K4me1")){
    return("1kb")
  }else{
    return("10kb")
  }
}


volcano_plot <- function(tissue){
  if(tissue == "ileum"){
    antibodys <- c("H3K27me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
  }else{
    antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
  }
  colour<- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p_list <- list(list(),list())
  i <- 1
  for (antibody in antibodys){
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff.csv"))
    p_list[[1]][[i]] <- ggplot(
      df, aes(x = LogFC.old.young, y = -log10(FDR.old.young))) +
      geom_point(aes(color = Significant), size=2,show.legend = FALSE) +
      scale_color_manual(values=colour) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
      ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (p-value)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 30))+
      annotate("text", x = min(df$LogFC.old.young), y = max(-log10(df$FDR.old.young)), label = nrow(df[which(df$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
      annotate("text", x = max(df$LogFC.old.young), y = max(-log10(df$FDR.old.young)), label = nrow(df[which(df$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
    if(antibody %in% c("H3K27ac","H3K4me3","H3K4me1")){
      peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size(antibody),"_in_young_old_macs_narrowpeak.bed"))  
    }else{
      peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size(antibody),"_in_young_old_merge-W1000-G3000-E100.bed"))
    }
    df <- df[which(df$Geneid %in% peak$V4),]
    p_list[[2]][[i]] <- ggplot(
      df, aes(x = LogFC.old.young, y = -log10(FDR.old.young))) +
      geom_point(aes(color = Significant), size=2,show.legend = FALSE) +
      scale_color_manual(values=colour) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
      ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (p-value)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 30))+
      annotate("text", x = min(df$LogFC.old.young), y = max(-log10(df$FDR.old.young)), label = nrow(df[which(df$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
      annotate("text", x = max(df$LogFC.old.young), y = max(-log10(df$FDR.old.young)), label = nrow(df[which(df$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
    i <- i+1
    }
  combined_plot <- plot_a_list(p_list[[1]],no_of_rows = 2,no_of_cols = 3)
  ggsave(paste0("result/",tissue,"/all_diff_volcano_plot.png"),combined_plot,width = 18,height = 12,type="cairo")
  combined_plot <- plot_a_list(p_list[[2]],no_of_rows = 2,no_of_cols = 3)
  ggsave(paste0("result/",tissue,"/all_diff_volcano_plot_in_peak.png"),combined_plot,width = 18,height = 12,type="cairo")
}
tissue <- "thymus"
volcano_plot(tissue)
  