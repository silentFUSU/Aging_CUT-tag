rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(edgeR)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT"))
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

diff_analysis_batch_or_not_bothsig_onlysig <- function(tissue,antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  if (!file.exists(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))) {  
    return(0)  
  }  
  batch <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  not_batch <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff.csv"))
  df <- merge(batch[,c("Geneid","Chr","Start","End","FDR.old.young")],not_batch[,c("Geneid","FDR.old.young")],by="Geneid")
  colnames(df)[c(5,6)]<-c("batch","not_batch")
  df$significant_condition <- "both significant"
  df$significant_condition[which(df$batch < 0.05 & df$not_batch > 0.05)] <-"only ~age+batch significant"
  df$significant_condition[which(df$batch > 0.05 & df$not_batch < 0.05)] <- "only ~age significant"
  df$significant_condition[which(df$batch > 0.05 & df$not_batch > 0.05)] <- "insignificant"
  df <- df[which(df$significant_condition != "insignificant"),]
  to_plot <- as.data.frame(table(df$significant_condition))
  to_plot$Percentage <- to_plot$Freq/sum(to_plot$Freq)*100
  # to_plot$Label <- paste0(to_plot$Var1, " (", round(to_plot$Percentage, 1), "%)")  
  colors <- read.table("data/samples/7_distinct_color.txt")
  colors <-setNames(colors$V1,c("insignificant","both significant","only ~age+batch significant","only ~age significant"))
  p <- ggplot(to_plot, aes(x = "", y = Freq, fill = Var1)) +  
    geom_bar(width = 1, stat = "identity", color = "white") +  
    scale_fill_manual(values = colors)+
    coord_polar("y", start = 0) +  
    theme_void() + 
    labs(fill = NULL) +  
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody)) + 
    theme(legend.position = "right",plot.title = element_text(hjust = 0.5),text = element_text(size = 16))
  dir.create(paste0("data/samples/all/batch_remove_or_keep/",antibody),showWarnings = F)
  write.csv(df,paste0("data/samples/all/batch_remove_or_keep/",antibody,"/",tissue,"_batch_remove_or_keep_significant.csv"))
  return(p)
  }
p_list <- list()
antibody <- "H3K9me3"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  p_list[[i]] <- diff_analysis_batch_or_not_bothsig_onlysig(tissue,antibody)
}
p_list <- p_list[-19]
combined_plot <- plot_a_list(p_list, 4, 7)
ggsave("result/all/all_tissues_batch_or_not_pie.png",combined_plot,width = 35,height = 20,type="cairo")
