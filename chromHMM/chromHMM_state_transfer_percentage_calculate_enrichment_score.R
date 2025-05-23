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

transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))

background_transfer_matrix <- data.frame(tissue = character(),
                              young_state = character(),
                              old_state = character(),
                              Freq = numeric(),
                              stringsAsFactors = FALSE)

tissues <-   c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
               "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")
state_num <- "11"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste0(young1$V1,"-",young1$V2,"-",young1$V3)
  young2$label <- paste0(young2$V1,"-",young2$V2,"-",young2$V3)
  states <- sort(unique(young1$V4))
  for(state in states){
    state_region <- young1[which(young1$V4 == state),]
    transfer_state <- young2[which(young2$label %in% state_region$label),]
    t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
    colnames(t_transfer_matrix) <- c("young2_state","Freq")
    t_transfer_matrix <- data.frame(tissue = tissue_label_change(tissue), 
                                    young1_state = state, 
                                    t_transfer_matrix)
    background_transfer_matrix <- rbind(background_transfer_matrix,t_transfer_matrix)
  }
}
transfer_matrix$young_state <- factor(transfer_matrix$young_state,levels=c(paste0("E",c(1:length(states)))))
transfer_matrix$old_state <- factor(transfer_matrix$old_state,levels=c(paste0("E",c(1:length(states)))))

background_transfer_matrix$young1_state <- factor(background_transfer_matrix$young1_state, levels=c(paste0("E",c(1:length(states)))))
background_transfer_matrix$young2_state <- factor(background_transfer_matrix$young2_state, levels=c(paste0("E",c(1:length(states)))))

result_transfer_matrix <- transfer_matrix %>%
  group_by(young_state, old_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_transfer_matrix2 <- result_transfer_matrix 
result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 

result_background_transfer_matrix <- background_transfer_matrix %>%
  group_by(young1_state, young2_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_background_transfer_matrix2 <- result_background_transfer_matrix
result_background_transfer_matrix2$percent_freq <- result_background_transfer_matrix2$mean_freq/sum(result_background_transfer_matrix2$mean_freq)*100
result_background_transfer_matrix_to_plot <- dcast(result_background_transfer_matrix2, formula = young1_state~young2_state, value.var = "percent_freq") 

rownames(result_transfer_matrix_to_plot) <- result_transfer_matrix_to_plot$young_state
result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]

rownames(result_background_transfer_matrix_to_plot) <- result_background_transfer_matrix_to_plot$young1_state
result_background_transfer_matrix_to_plot <- result_background_transfer_matrix_to_plot[,-1]

result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
labels <-as.character(rownames(result_to_plot))

pheatmap::pheatmap(result_to_plot,
                   cluster_rows = F,cluster_cols = F, 
                   breaks = seq(0, 5, length.out = 101),
                   display_numbers = T,labels_row = labels,
                   labels_col = labels,fontsize = 10,na_col = "white")

pheatmap::pheatmap(result_background_transfer_matrix_to_plot,
                   cluster_rows = F,cluster_cols = F, 
                   breaks = seq(0, 5, length.out = 101),
                   display_numbers = T,labels_row = labels,
                   labels_col = labels,fontsize = 10,na_col = "white")



###### split by tissues
tissues <-   c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
               "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")
state_num <- "11"
for(tissue in tissues){
  tissue_label <- tissue_label_change(tissue)
  transfer <- transfer_matrix[which(transfer_matrix$tissue==tissue_label),]
  background <- background_transfer_matrix[which(background_transfer_matrix$tissue==tissue_label),]
  result_transfer_matrix <- transfer  %>%
    group_by(young_state, old_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_transfer_matrix2 <- result_transfer_matrix 
  result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
  result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
  result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
  
  result_background_transfer_matrix <- background %>%
    group_by(young1_state, young2_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_background_transfer_matrix2 <- result_background_transfer_matrix
  result_background_transfer_matrix2$percent_freq <- result_background_transfer_matrix2$mean_freq/sum(result_background_transfer_matrix2$mean_freq)*100
  result_background_transfer_matrix_to_plot <- dcast(result_background_transfer_matrix2, formula = young1_state~young2_state, value.var = "percent_freq") 
  
  rownames(result_transfer_matrix_to_plot) <- result_transfer_matrix_to_plot$young_state
  result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
  
  rownames(result_background_transfer_matrix_to_plot) <- result_background_transfer_matrix_to_plot$young1_state
  result_background_transfer_matrix_to_plot <- result_background_transfer_matrix_to_plot[,-1]
  
  result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
  labels <-as.character(rownames(result_to_plot))
  result_to_plot[is.na(result_to_plot)] <- 0
  pheatmap::pheatmap(result_to_plot,
                     cluster_rows = F,cluster_cols = F, 
                     breaks = seq(0, 5, length.out = 101),
                     display_numbers = T,labels_row = labels,
                     labels_col = labels,fontsize = 10,na_col = "white",main = tissue_label,
                     filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/",tissue,"_state_transition_enrichment_score_young_replicates.png"),width = 5,height = 5)
  }

#### use enrichment score clustering tissues
summary <- data.frame()
for(tissue in tissues){
  tissue_label <- tissue_label_change(tissue)
  transfer <- transfer_matrix[which(transfer_matrix$tissue==tissue_label),]
  background <- background_transfer_matrix[which(background_transfer_matrix$tissue==tissue_label),]
  result_transfer_matrix <- transfer  %>%
    group_by(young_state, old_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_transfer_matrix2 <- result_transfer_matrix 
  result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
  result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
  result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
  
  result_background_transfer_matrix <- background %>%
    group_by(young1_state, young2_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_background_transfer_matrix2 <- result_background_transfer_matrix
  result_background_transfer_matrix2$percent_freq <- result_background_transfer_matrix2$mean_freq/sum(result_background_transfer_matrix2$mean_freq)*100
  result_background_transfer_matrix_to_plot <- dcast(result_background_transfer_matrix2, formula = young1_state~young2_state, value.var = "percent_freq") 
  
  rownames(result_transfer_matrix_to_plot) <- result_transfer_matrix_to_plot$young_state
  result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
  
  rownames(result_background_transfer_matrix_to_plot) <- result_background_transfer_matrix_to_plot$young1_state
  result_background_transfer_matrix_to_plot <- result_background_transfer_matrix_to_plot[,-1]
  
  result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
  labels <-as.character(rownames(result_to_plot))
  result_to_plot[is.na(result_to_plot)] <- 0
  result_to_plot$young_state <- rownames(result_to_plot)
  result <- melt(result_to_plot,value.name = "proportion")
  result$young_state <- factor(result$young_state,levels = c(paste0("E",c(1:length(states)))))
  result$variable <- factor(result$variable,levels = c(paste0("E",c(1:length(states)))))
  result <- result %>%
    arrange(young_state, variable)
  result$label <- paste0(result$young_state,"-",result$variable)
  result <- result[,c("label","proportion")]
  colnames(result)[2] <- tissue_label_change(tissue)
  if(nrow(summary)==0){
    summary <- result
  }else{
    summary <- merge(summary,result,by="label")
  }
}
combinations <- expand.grid(paste0("E",c(1:length(states))),paste0("E",c(1:length(states))))
combinations$Var1 <- factor(combinations$Var1,paste0("E",c(1:length(states))))
combinations$Var2 <- factor(combinations$Var2,paste0("E",c(1:length(states))))
combinations<- combinations %>%
  arrange(Var1, Var2)
combinations <- apply(combinations, 1, function(x) paste(x, collapse = "-"))
summary$label <- factor(summary$label,levels = combinations)
summary <- summary[order(summary$label),]
rownames(summary) <- summary$label
summary <- summary[,-1]
pheatmap::pheatmap(summary,
                   breaks = seq(0, 5, length.out = 101),
                   fontsize = 10,na_col = "white",
                   filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/all_tissue_state_transition_enrichment_score_young_replicates_heatmap.png"),
                   width = 10,height = 20
                   )

cor_to_plot <- cor(summary,method = "spearman")
pheatmap::pheatmap(cor_to_plot,breaks = seq(-1, 1, length.out = 101),
                   filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/all_tissue_state_transition_enrichment_score_young_replicates_correaltion.png"),
                   width = 5,height = 5)
