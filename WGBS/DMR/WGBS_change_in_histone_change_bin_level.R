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
                        nrow = no_of_rows, ncol = no_of_cols)
}
# tissues <- c("liver","lung","kidney","ileum","Hip","mammarygland","skin","bonemarrow","jejunum","colon","ovary","CB")
# tissues <- c("BAT","thymus","testis")
tissues <- sort(c("liver","lung","kidney","ileum","Hip","mammarygland","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis","heart","muscle","stomach","bladder","aorta","tongue","spleen","pancreas","brain","cecum","uterus","iWAT"))

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

WGBS_change_in_histone_change <- function(tissue,antibody){
  diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff.csv"))
  bin_in_peaks <- read.table(paste0("data/samples/",tissue,"/H3K36me3/bed/H3K36me3_10kb_in_young_old_merge-W1000-G3000-E100.bed"))
  diff <- diff[which(diff$Geneid %in% bin_in_peaks$V4),]
  increase <- diff$Geneid[which(diff$Significant=="Up")]
  decrease <- diff$Geneid[which(diff$Significant=="Down")]
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size(antibody),"_bins_all_depth.csv"))
  df <- df[which(df$total_V5>15),]
  df <- merge(df,search_table[,c(3:5)],by = "sample")
  df$percent <- df$percent*100
  df$age[which(df$age=="3M")] <- "young"
  df$age[which(df$age=="24M")] <- "old"
  df$age <- factor(df$age,levels=c("young","old"))
  increase_bin <- df[which(df$label %in% increase),]
  decrease_bin <- df[which(df$label %in% decrease),]
  p_list <- list()
  if(nrow(increase_bin) > 0){
    t <- t.test(increase_bin$percent[which(increase_bin$age=="old")],increase_bin$percent[which(increase_bin$age=="young")])
    increase_plot <- ggplot(increase_bin, aes(x = age, y = percent,fill=sample)) +  
      scale_fill_brewer(palette = "Pastel1") +
      geom_boxplot() +  
      theme_minimal() +
      theme(text = element_text(size = 20),legend.position = "none") +
      labs(title = paste0(tissue_label_change(tissue),"\n",antibody," Increase Region"), x = NULL, y = "CG%") +
      annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
               hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
      annotate("text", x = -Inf, y = -Inf, label = paste("bin number =",  length(unique(increase_bin$label))  ),   
               hjust = 0, vjust = -1.1, size = 5, colour = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
               hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
      annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
               hjust = 0, vjust = 1.1, size = 5, colour = "red")
    p_list[["increase"]] <- increase_plot
  }
  if(nrow(decrease_bin) > 0){
    t <- t.test(decrease_bin$percent[which(decrease_bin$age=="old")],decrease_bin$percent[which(decrease_bin$age=="young")])
    decrease_plot <- ggplot(decrease_bin, aes(x = age, y = percent,fill=sample)) +  
      scale_fill_brewer(palette = "Pastel1") +
      geom_boxplot() +  
      theme_minimal() +
      theme(text = element_text(size = 20),legend.position = "none") +
      labs(title = paste0(tissue_label_change(tissue),"\n",antibody," Decrease Region"), x = NULL, y = "CG%") +
      annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
               hjust = 1.1, vjust = -1.1, size = 5, colour = "red") +
      annotate("text", x = -Inf, y = -Inf, label = paste("bin number =",  length(unique(decrease_bin$label))  ),   
               hjust = 0, vjust = -1.1, size = 5, colour = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
               hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
      annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
               hjust = 0, vjust = 1.1, size = 5, colour = "red")
    p_list[["decrease"]] <- decrease_plot
  }
  return(p_list)
}
  
p_list <- list()
i <- 1
for(tissue in tissues){
  t_p_list <- WGBS_change_in_histone_change(tissue,"H3K4me1")
  if ("increase" %in% names(t_p_list)) {
      p_list[[i]] <- t_p_list[["increase"]]
      i <- i+1
  }
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 7)
ggsave("result/WGBS/all_tissues_change_in_H3K4me1_increase.png",combined_plot,width = 35,height = 28,type="cairo")
