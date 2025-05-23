rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(reshape2)
library(ggalluvial)  
library(bubbleHeatmap)
options(bitmapType = "cairo")  
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
      tissue_label <- "Mammarygland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
state_num <- "11"
#### use enrichment score clustering tissues
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_cross_comparison.csv"),row.names = 1)
background_transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_backgroud_cross_comparison.csv"),row.names = 1)
transfer_matrix$young_state <- factor(transfer_matrix$young_state,levels=c(paste0("E",c(1:state_num))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state,levels=c(paste0("E",c(1:state_num))))
background_transfer_matrix$rep1_state <- factor(background_transfer_matrix$rep1_state, levels=c(paste0("E",c(1:state_num))))
background_transfer_matrix$rep2_state <- factor(background_transfer_matrix$rep2_state, levels=c(paste0("E",c(1:state_num))))

result_transfer_matrix <- transfer_matrix %>%
  group_by(young_state, old_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_transfer_matrix2 <- result_transfer_matrix 
result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 

result_background_transfer_matrix <- background_transfer_matrix %>%
  group_by(rep1_state, rep2_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_background_transfer_matrix2 <- result_background_transfer_matrix
result_background_transfer_matrix2$percent_freq <- result_background_transfer_matrix2$mean_freq/sum(result_background_transfer_matrix2$mean_freq)*100
result_background_transfer_matrix_to_plot <- dcast(result_background_transfer_matrix2, formula = rep1_state~rep2_state, value.var = "percent_freq") 

rownames(result_transfer_matrix_to_plot) <- result_transfer_matrix_to_plot$young_state
result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]

rownames(result_background_transfer_matrix_to_plot) <- result_background_transfer_matrix_to_plot$young1_state
result_background_transfer_matrix_to_plot <- result_background_transfer_matrix_to_plot[,-1]

result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
labels <-as.character(rownames(result_to_plot))

diag(result_to_plot) <- NA
display_numbers_matrix <- matrix(nrow = nrow(result_to_plot), ncol = ncol(result_to_plot))
for (i in 1:nrow(result_to_plot)) {
  for (j in 1:ncol(result_to_plot)) {
    if (is.na(result_to_plot[i, j])) {
      display_numbers_matrix[i, j] <- ""
    } else {
      display_numbers_matrix[i, j] <- format(result_to_plot[i, j], nsmall = 2, digits = 2)
    }
  }
}
pheatmap::pheatmap(result_to_plot,
                   cluster_rows = F,cluster_cols = F, 
                   breaks = seq(0, 5, length.out = 101),
                   display_numbers = display_numbers_matrix,
                   labels_row = labels,
                   labels_col = labels,fontsize = 10,na_col = "white")

result_transfer_matrix2 <- result_transfer_matrix 
result_transfer_matrix2$mean_freq[which(result_transfer_matrix2$young_state==result_transfer_matrix2$old_state)] <- NA
result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq,na.rm = T)*100
bubble_result_transfer_matrix_to_plot <- reshape2::dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
rownames(bubble_result_transfer_matrix_to_plot) <- bubble_result_transfer_matrix_to_plot$young_state
bubble_result_transfer_matrix_to_plot <- bubble_result_transfer_matrix_to_plot[,-1]

# result_to_plot[is.na(result_to_plot)] <- 0
# bubble_result_transfer_matrix_to_plot[is.na(bubble_result_transfer_matrix_to_plot)] <- 0
tree <-  bubbleHeatmap(as.matrix(result_to_plot), as.matrix(bubble_result_transfer_matrix_to_plot),
                       colorLim = c(0,4),sizeLim = c(0,10),
                       legendTitles = c("Proportion", "Enrichment score"))
grid.newpage()
grid.draw(tree)
