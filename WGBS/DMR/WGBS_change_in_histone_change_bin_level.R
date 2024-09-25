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
library(DSS)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
tissues <- c("liver","lung","kidney","ileum","Hip")
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
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return("10kb")
  }else{
    return("1kb")
  }
}

search_table <-read.csv("data/samples/all/WGBS_search_table.csv")
colnames(search_table)[3] <- "sample"
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")


for(tissue in tissues){
  p_list <- list()
  i=1
  if(tissue == "ileum"){
    antibodys <- c("H3K27me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
  }else{
    antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
  }
  for(antibody in antibodys){
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff.csv"))
    increase <- diff$Geneid[which(diff$Significant=="Up")]
    decrease <- diff$Geneid[which(diff$Significant=="Down")]
    df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size(antibody),"_bins_all_depth.csv"))
    df <- df[which(df$total_V5>15),]
    df <- merge(df,search_table[,c(3:5)],by = "sample")
    df$percent <- df$percent*100
    df$age <- factor(df$age,levels=c("3M","24M"))
    increase_bin <- df[which(df$label %in% increase),]
    decrease_bin <- df[which(df$label %in% decrease),]
    if(nrow(increase_bin) > 0){
      t <- t.test(increase_bin$percent[which(increase_bin$age=="3M")],increase_bin$percent[which(increase_bin$age=="24M")])
      increase_plot <- ggplot(increase_bin, aes(x = age, y = percent,fill=sample)) +  
        scale_fill_brewer(palette = "Pastel1") +
        geom_boxplot() +  
        theme_minimal() +
        theme(text = element_text(size = 20)) +
        labs(title = paste0(tissue_label_change(tissue),"\n",antibody," Increase Region"), x = NULL, y = "CG%") +
        annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
                 hjust = 1.1, vjust = -1.1, size = 5, colour = "red")
      p_list[[i]] <- increase_plot
      i <- i+1
    }
    if(nrow(decrease_bin) > 0){
      t <- t.test(decrease_bin$percent[which(decrease_bin$age=="3M")],decrease_bin$percent[which(decrease_bin$age=="24M")])
      decrease_plot <- ggplot(decrease_bin, aes(x = age, y = percent,fill=sample)) +  
        scale_fill_brewer(palette = "Pastel1") +
        geom_boxplot() +  
        theme_minimal() +
        theme(text = element_text(size = 20)) +
        labs(title = paste0(tissue_label_change(tissue),"\n",antibody," Decrease Region"), x = NULL, y = "CG%") +
        annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
                 hjust = 1.1, vjust = -1.1, size = 5, colour = "red")
      p_list[[i]] <- decrease_plot
      i <- i+1
    }
  }
  combined_plot <- plot_a_list(p_list,no_of_cols = 2, no_of_rows = length(p_list)/2)
  ggsave(paste0("result/WGBS/",tissue,"/WGBS_change_in_histone_change/plot/all_antibodys_box_plot_without_removing_batch_effect.png"),width = 12,height = 6*length(p_list)/2,type="cairo")
}
