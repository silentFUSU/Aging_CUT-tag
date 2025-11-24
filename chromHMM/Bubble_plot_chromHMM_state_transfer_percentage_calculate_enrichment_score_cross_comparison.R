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
state_num <- 15
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
#### use enrichment score clustering tissues  
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/state_transfer/state_transfer_cross_comparison.csv"),row.names = 1)
transfer_matrix$young_state <- values[match(transfer_matrix$young_state, keys)]
transfer_matrix$old_state <- values[match(transfer_matrix$old_state, keys)]
transfer_matrix$young_state <- paste0("E",transfer_matrix$young_state)
transfer_matrix$old_state <- paste0("E",transfer_matrix$old_state)


background_transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/state_transfer/state_transfer_backgroud_cross_comparison.csv"),row.names = 1)
background_transfer_matrix$rep1_state <- values[match(background_transfer_matrix$rep1_state, keys)]
background_transfer_matrix$rep2_state <- values[match(background_transfer_matrix$rep2_state, keys)]
background_transfer_matrix$rep1_state <- paste0("E",background_transfer_matrix$rep1_state)
background_transfer_matrix$rep2_state <- paste0("E",background_transfer_matrix$rep2_state)


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

result_to_plot <- log2(result_to_plot)
result_transfer_matrix2 <- result_transfer_matrix 
result_transfer_matrix2$mean_freq[which(result_transfer_matrix2$young_state==result_transfer_matrix2$old_state)] <- NA
result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq,na.rm = T)*100
bubble_result_transfer_matrix_to_plot <- reshape2::dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
rownames(bubble_result_transfer_matrix_to_plot) <- bubble_result_transfer_matrix_to_plot$young_state
bubble_result_transfer_matrix_to_plot <- bubble_result_transfer_matrix_to_plot[,-1]

# result_to_plot[is.na(result_to_plot)] <- 0
# bubble_result_transfer_matrix_to_plot[is.na(bubble_result_transfer_matrix_to_plot)] <- 0
customColors <- colorRampPalette(c("blue", "white", "red"))
tree <-  bubbleHeatmap(as.matrix(result_to_plot), as.matrix(bubble_result_transfer_matrix_to_plot),
                       colorLim = c(-1,1),sizeLim = c(0,10),
                       legendTitles = c("Proportion", "log2(Enrichment score)"),colorSeq = customColors(100),xTitle = "Old state",yTitle = "Young state")

svglite::svglite("result/figures/chromHMM_15_bubble_plot.svg", width = 8, height = 6)

grid.newpage()
grid.draw(tree)
dev.off()

### each tissue
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))

for(tissue in tissues){
  state_num <- 15
  transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_cross_comparison.csv"),row.names = 1)
  background_transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_backgroud_cross_comparison.csv"),row.names = 1)
  transfer_matrix$young_state <- factor(transfer_matrix$young_state,levels=c(paste0("E",c(1:state_num))))
  transfer_matrix$old_state <- factor(transfer_matrix$old_state,levels=c(paste0("E",c(1:state_num))))
  background_transfer_matrix$rep1_state <- factor(background_transfer_matrix$rep1_state, levels=c(paste0("E",c(1:state_num))))
  background_transfer_matrix$rep2_state <- factor(background_transfer_matrix$rep2_state, levels=c(paste0("E",c(1:state_num))))
  
  transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue == tissue_label_change(tissue)),]
  background_transfer_matrix <- background_transfer_matrix[which(background_transfer_matrix$tissue == tissue_label_change(tissue)),]
  
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
  
  result_transfer_matrix_to_plot[is.na(result_transfer_matrix_to_plot)] <- 0
  result_background_transfer_matrix_to_plot[is.na(result_background_transfer_matrix_to_plot)] <- 0
  
  result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
  labels <-as.character(rownames(result_to_plot))
  
  diag(result_to_plot) <- NA

  result_to_plot <- as.matrix(result_to_plot)
  finite_values <- result_to_plot
  finite_values[is.infinite(finite_values) | is.na(finite_values)] <- NA
  max_value <- max(finite_values, na.rm = TRUE)
  result_to_plot[is.infinite(result_to_plot)] <- max_value
  result_to_plot <- as.data.frame(result_to_plot)
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
  # pheatmap::pheatmap(result_to_plot,
  #                    cluster_rows = F,cluster_cols = F, 
  #                    breaks = seq(0, 5, length.out = 101),
  #                    display_numbers = display_numbers_matrix,
  #                    labels_row = labels,
  #                    labels_col = labels,fontsize = 10,na_col = "white")
  
  result_transfer_matrix2 <- result_transfer_matrix 
  result_transfer_matrix2$mean_freq[which(result_transfer_matrix2$young_state==result_transfer_matrix2$old_state)] <- NA
  result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq,na.rm = T)*100
  bubble_result_transfer_matrix_to_plot <- reshape2::dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
  rownames(bubble_result_transfer_matrix_to_plot) <- bubble_result_transfer_matrix_to_plot$young_state
  bubble_result_transfer_matrix_to_plot <- bubble_result_transfer_matrix_to_plot[,-1]
  
  customColors <- colorRampPalette(c("#4589C8FF", "white", "#EE7C7AFF"))
  tree <-  bubbleHeatmap(as.matrix(result_to_plot), as.matrix(bubble_result_transfer_matrix_to_plot),
                         colorLim = c(0,2),sizeLim = c(0,10),
                         legendTitles = c("Proportion", "Enrichment score"),colorSeq = customColors(100),xTitle = "Old state",yTitle = "Young state",plotTitle = tissue_label_change(tissue))
  
  svglite::svglite(paste0("result/Sup_figures/chromHMM_15_each_tissue_bubble/",tissue,"_chromHMM_15_bubble_plot.svg"), width = 8, height = 6)
  grid.newpage()
  grid.draw(tree)
  dev.off()
  
}
