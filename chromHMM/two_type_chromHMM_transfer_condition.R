rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(corrplot)
options(scipen = 999)  
state_num <- "14"
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))


type1 <- c("Liver","Bone Marrow","Heart","Skin","Spleen","Cecum","Colon","Lung","Cortex","Hippocampus","Aorta","Muscle","Stomach")
type2 <- c("Ovary","Mammarygland","Tongue","Uterus","Thymus","Jejunum","Testis","iWAT","BAT","Kidney","Cerebellum","Bladder","Pancreas")
# scale by row
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue %in% type2),]
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:state_num))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:state_num))))


result_transfer_matrix <- transfer_matrix %>%
  group_by(young_state, old_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_transfer_matrix2 <- result_transfer_matrix %>%
  group_by(young_state) %>%
  mutate(percent_freq = mean_freq/sum(mean_freq))
result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
labels <-as.character(result_transfer_matrix_to_plot$young_state)
result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
result_transfer_matrix_to_plot[is.na(result_transfer_matrix_to_plot)] <- 0
pheatmap::pheatmap(result_transfer_matrix_to_plot,
                   cluster_rows = F,cluster_cols = F, 
                   breaks = seq(0, 0.2, length.out = 101),
                   display_numbers = T,labels_row = labels,
                   labels_col = labels,fontsize = 10)

#Scale by all condition, remove empty states
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix <- transfer_matrix[-which(transfer_matrix$young_state %in% c("E3","E12") | transfer_matrix$old_state %in% c("E3","E12")),]
transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue %in% type2),]
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:state_num))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:state_num))))

result_transfer_matrix <- transfer_matrix %>%
  group_by(young_state, old_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_transfer_matrix2 <- result_transfer_matrix 
result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100

result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
labels <-as.character(result_transfer_matrix_to_plot$young_state)
result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
result_transfer_matrix_to_plot[is.na(result_transfer_matrix_to_plot)] <- 0
pheatmap::pheatmap(result_transfer_matrix_to_plot,
                   cluster_rows = F,cluster_cols = F, 
                   breaks = seq(0, 5, length.out = 101),
                   display_numbers = T,labels_row = labels,
                   labels_col = labels,fontsize = 10,na_col = "white")

