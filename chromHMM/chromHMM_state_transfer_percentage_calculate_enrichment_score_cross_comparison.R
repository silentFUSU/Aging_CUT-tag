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
state_num <- 15
background_transfer_matrix <- data.frame(tissue = character(),
                                         rep1_state = character(),
                                         rep2_state = character(),
                                         Freq = numeric(),
                                         stringsAsFactors = FALSE)
tissues <-   c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
               "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste0(young1$V1,"-",young1$V2,"-",young1$V3)
  young2$label <- paste0(young2$V1,"-",young2$V2,"-",young2$V3)
  states <- sort(unique(young1$V4))
  young_background_transfer_matrix <- data.frame()
  for(state in states){
    state_region <- young1[which(young1$V4 == state),]
    transfer_state <- young2[which(young2$label %in% state_region$label),]
    t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
    colnames(t_transfer_matrix) <- c("young2_state","Freq")
    t_transfer_matrix <- data.frame(tissue = tissue_label_change(tissue), 
                                    young1_state = state, 
                                    t_transfer_matrix)
    young_background_transfer_matrix <- rbind(young_background_transfer_matrix,t_transfer_matrix)
  }
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste0(old1$V1,"-",old1$V2,"-",old1$V3)
  old2$label <- paste0(old2$V1,"-",old2$V2,"-",old2$V3)
  states <- sort(unique(old1$V4))
  old_background_transfer_matrix <- data.frame()
  for(state in states){
    state_region <- old1[which(old1$V4 == state),]
    transfer_state <- old2[which(old2$label %in% state_region$label),]
    t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
    colnames(t_transfer_matrix) <- c("old2_state","Freq")
    t_transfer_matrix <- data.frame(tissue = tissue_label_change(tissue), 
                                    old1_state = state, 
                                    t_transfer_matrix)
    old_background_transfer_matrix <- rbind(old_background_transfer_matrix,t_transfer_matrix)
  }
  young_background_transfer_matrix$label <- paste0(young_background_transfer_matrix$young1_state,"-",young_background_transfer_matrix$young2_state)
  old_background_transfer_matrix$label <- paste0(old_background_transfer_matrix$old1_state,"-",old_background_transfer_matrix$old2_state)
  t_background_transfer_matrix <- merge(young_background_transfer_matrix, old_background_transfer_matrix, by="label", all=T)
  t_background_transfer_matrix <- data.frame(tissue=tissue_label_change(tissue),
                                             label=t_background_transfer_matrix$label,
                                             Freq= rowMeans(t_background_transfer_matrix[,c("Freq.x","Freq.y")], na.rm = T))
  t_background_transfer_matrix <- t_background_transfer_matrix %>%
    separate(
      col = label,      
      into = c("rep1_state", "rep2_state"),  
      sep = "-",            
      remove = TRUE,      
      convert = FALSE
    )
  background_transfer_matrix <- rbind(background_transfer_matrix,t_background_transfer_matrix)
}
dir.create(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/state_transfer/"))
write.csv(background_transfer_matrix,paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/state_transfer/state_transfer_backgroud_cross_comparison.csv"))

transfer_matrix <- data.frame(tissue = character(),
                              young_state = character(),
                              old_state = character(),
                              Freq = numeric(),
                              stringsAsFactors = FALSE)
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste0(young1$V1,"-",young1$V2,"-",young1$V3)
  young2$label <- paste0(young2$V1,"-",young2$V2,"-",young2$V3)
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste0(old1$V1,"-",old1$V2,"-",old1$V3)
  old2$label <- paste0(old2$V1,"-",old2$V2,"-",old2$V3)
  young <- list(young1,young2)
  old <- list(old1,old2)
  tissue_transfer_matrix <- list()
  for(j in c(1:2)){
    t_young <- young[[j]]
    t_old <- old[[j]]
    tissue_transfer_matrix[[j]] <- data.frame()
    t_young <- t_young[which(t_young$label %in% t_old$label),]
    t_old <-  t_old[which(t_old$label %in% t_young$label),]
    states <- sort(unique(t_young$V4))
    for(state in states){
      state_region <- t_young[which(t_young$V4 == state),]
      transfer_state <- t_old[which(t_old$label %in% state_region$label),]
      t_tissue_transfer_matrix <- as.data.frame(table(transfer_state$V4))
      colnames(t_tissue_transfer_matrix) <- c("old_state","Freq")
      t_tissue_transfer_matrix <- data.frame(tissue = tissue_label_change(tissue),
                                             young_state = state,
                                             t_tissue_transfer_matrix)
      t_tissue_transfer_matrix$label <- paste0(t_tissue_transfer_matrix$young_state,"-",t_tissue_transfer_matrix$old_state)
      tissue_transfer_matrix[[j]] <- rbind(tissue_transfer_matrix[[j]],t_tissue_transfer_matrix)
      }
    }
  extract_columns <- lapply(tissue_transfer_matrix, function(df) {
    select(df, label, Freq)
  })
  tissue_transfer_matrix_summary <- Reduce(function(x, y) {
    full_join(x, y, by = "label", suffix = c("", ".y"))
  }, extract_columns)
  
  tissue_transfer_matrix_summary$average_Freq <- rowMeans(tissue_transfer_matrix_summary[,-1],na.rm = T)
  tissue_transfer_matrix_summary <- tissue_transfer_matrix_summary[,c("label","average_Freq")]
  tissue_transfer_matrix_summary <- tissue_transfer_matrix_summary %>%
    separate(
      col = label,      
      into = c("young_state", "old_state"),  
      sep = "-",            
      remove = TRUE,      
      convert = FALSE
    )
  tissue_transfer_matrix_summary$tissue <- tissue_label_change(tissue)
  colnames(tissue_transfer_matrix_summary)[3] <- "Freq"
  tissue_transfer_matrix_summary <- tissue_transfer_matrix_summary[,c("tissue","young_state","old_state","Freq")]
  transfer_matrix <- rbind(transfer_matrix, tissue_transfer_matrix_summary)
}
write.csv(transfer_matrix,paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/state_transfer/state_transfer_cross_comparison.csv"))

# transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_cross_comparison.csv"),row.names = 1)
# background_transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_backgroud_cross_comparison.csv"),row.names = 1)

# transfer_matrix$young_state <- factor(transfer_matrix$young_state,levels=c(paste0("E",c(1:state_num))))
# transfer_matrix$old_state <- factor(transfer_matrix$old_state,levels=c(paste0("E",c(1:state_num))))
# background_transfer_matrix$rep1_state <- factor(background_transfer_matrix$rep1_state, levels=c(paste0("E",c(1:state_num))))
# background_transfer_matrix$rep2_state <- factor(background_transfer_matrix$rep2_state, levels=c(paste0("E",c(1:state_num))))
# 
# result_transfer_matrix <- transfer_matrix %>%
#   group_by(young_state, old_state) %>%
#   summarise(mean_freq = mean(Freq, na.rm = TRUE))
# result_transfer_matrix2 <- result_transfer_matrix 
# result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
# result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
# 
# result_background_transfer_matrix <- background_transfer_matrix %>%
#   group_by(rep1_state, rep2_state) %>%
#   summarise(mean_freq = mean(Freq, na.rm = TRUE))
# result_background_transfer_matrix2 <- result_background_transfer_matrix
# result_background_transfer_matrix2$percent_freq <- result_background_transfer_matrix2$mean_freq/sum(result_background_transfer_matrix2$mean_freq)*100
# result_background_transfer_matrix_to_plot <- dcast(result_background_transfer_matrix2, formula = rep1_state~rep2_state, value.var = "percent_freq") 
# 
# rownames(result_transfer_matrix_to_plot) <- result_transfer_matrix_to_plot$young_state
# result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
# 
# rownames(result_background_transfer_matrix_to_plot) <- result_background_transfer_matrix_to_plot$young1_state
# result_background_transfer_matrix_to_plot <- result_background_transfer_matrix_to_plot[,-1]
# 
# result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
# labels <-as.character(rownames(result_to_plot))
# 
# pheatmap::pheatmap(result_to_plot,
#                    cluster_rows = F,cluster_cols = F, 
#                    breaks = seq(0, 5, length.out = 101),
#                    display_numbers = T,labels_row = labels,
#                    labels_col = labels,fontsize = 10,na_col = "white")
# 
# pheatmap::pheatmap(result_background_transfer_matrix_to_plot,
#                    cluster_rows = F,cluster_cols = F, 
#                    breaks = seq(0, 5, length.out = 101),
#                    display_numbers = T,labels_row = labels,
#                    labels_col = labels,fontsize = 10,na_col = "white")
# 
# pheatmap::pheatmap(result_transfer_matrix_to_plot,
#                    cluster_rows = F,cluster_cols = F, 
#                    breaks = seq(0, 5, length.out = 101),
#                    display_numbers = T,labels_row = labels,
#                    labels_col = labels,fontsize = 10,na_col = "white")
# 
# ###### split by tissues
# tissues <-   c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
#                "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")
# state_num <- 11
# transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_cross_comparison.csv"),row.names = 1)
# background_transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_backgroud_cross_comparison.csv"),row.names = 1)
# transfer_matrix$young_state <- factor(transfer_matrix$young_state,levels=c(paste0("E",c(1:state_num))))
# transfer_matrix$old_state <- factor(transfer_matrix$old_state,levels=c(paste0("E",c(1:state_num))))
# background_transfer_matrix$rep1_state <- factor(background_transfer_matrix$rep1_state, levels=c(paste0("E",c(1:state_num))))
# background_transfer_matrix$rep2_state <- factor(background_transfer_matrix$rep2_state, levels=c(paste0("E",c(1:state_num))))
# 
# for(tissue in tissues){
#   tissue_label <- tissue_label_change(tissue)
#   transfer <- transfer_matrix[which(transfer_matrix$tissue==tissue_label),]
#   background <- background_transfer_matrix[which(background_transfer_matrix$tissue==tissue_label),]
#   result_transfer_matrix <- transfer  %>%
#     group_by(young_state, old_state) %>%
#     summarise(mean_freq = mean(Freq, na.rm = TRUE))
#   result_transfer_matrix2 <- result_transfer_matrix 
#   result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
#   result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
#   result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
#   
#   result_background_transfer_matrix <- background %>%
#     group_by(rep1_state, rep2_state) %>%
#     summarise(mean_freq = mean(Freq, na.rm = TRUE))
#   result_background_transfer_matrix2 <- result_background_transfer_matrix
#   result_background_transfer_matrix2$percent_freq <- result_background_transfer_matrix2$mean_freq/sum(result_background_transfer_matrix2$mean_freq)*100
#   result_background_transfer_matrix_to_plot <- dcast(result_background_transfer_matrix2, formula = rep1_state~rep2_state, value.var = "percent_freq") 
#   
#   rownames(result_transfer_matrix_to_plot) <- result_transfer_matrix_to_plot$young_state
#   result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
#   
#   rownames(result_background_transfer_matrix_to_plot) <- result_background_transfer_matrix_to_plot$young1_state
#   result_background_transfer_matrix_to_plot <- result_background_transfer_matrix_to_plot[,-1]
#   
#   result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
#   labels <-as.character(rownames(result_to_plot))
#   result_to_plot[is.na(result_to_plot)] <- 0
#   pheatmap::pheatmap(result_to_plot,
#                      cluster_rows = F,cluster_cols = F, 
#                      breaks = seq(0, 5, length.out = 101),
#                      display_numbers = T,labels_row = labels,
#                      labels_col = labels,fontsize = 10,na_col = "white",main = tissue_label,
#                      filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/",tissue,"_state_transition_enrichment_score_young_replicates_cross_comparison.png"),width = 5,height = 5)
# }
# 
# #### use enrichment score clustering tissues
# transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_cross_comparison.csv"),row.names = 1)
# background_transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer_backgroud_cross_comparison.csv"),row.names = 1)
# transfer_matrix$young_state <- factor(transfer_matrix$young_state,levels=c(paste0("E",c(1:state_num))))
# transfer_matrix$old_state <- factor(transfer_matrix$old_state,levels=c(paste0("E",c(1:state_num))))
# background_transfer_matrix$rep1_state <- factor(background_transfer_matrix$rep1_state, levels=c(paste0("E",c(1:state_num))))
# background_transfer_matrix$rep2_state <- factor(background_transfer_matrix$rep2_state, levels=c(paste0("E",c(1:state_num))))
# 
# summary <- data.frame()
# for(tissue in tissues){
#   tissue_label <- tissue_label_change(tissue)
#   transfer <- transfer_matrix[which(transfer_matrix$tissue==tissue_label),]
#   background <- background_transfer_matrix[which(background_transfer_matrix$tissue==tissue_label),]
#   result_transfer_matrix <- transfer  %>%
#     group_by(young_state, old_state) %>%
#     summarise(mean_freq = mean(Freq, na.rm = TRUE))
#   result_transfer_matrix2 <- result_transfer_matrix 
#   result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
#   result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
#   result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
#   
#   result_background_transfer_matrix <- background %>%
#     group_by(rep1_state, rep2_state) %>%
#     summarise(mean_freq = mean(Freq, na.rm = TRUE))
#   result_background_transfer_matrix2 <- result_background_transfer_matrix
#   result_background_transfer_matrix2$percent_freq <- result_background_transfer_matrix2$mean_freq/sum(result_background_transfer_matrix2$mean_freq)*100
#   result_background_transfer_matrix_to_plot <- dcast(result_background_transfer_matrix2, formula = rep1_state~rep2_state, value.var = "percent_freq") 
#   
#   rownames(result_transfer_matrix_to_plot) <- result_transfer_matrix_to_plot$young_state
#   result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
#   
#   rownames(result_background_transfer_matrix_to_plot) <- result_background_transfer_matrix_to_plot$rep1_state
#   result_background_transfer_matrix_to_plot <- result_background_transfer_matrix_to_plot[,-1]
#   
#   result_to_plot <- result_transfer_matrix_to_plot / result_background_transfer_matrix_to_plot
#   labels <-as.character(rownames(result_to_plot))
#   result_to_plot[is.na(result_to_plot)] <- 0
#   result_to_plot$young_state <- rownames(result_to_plot)
#   result <- melt(result_to_plot,value.name = "proportion")
#   result$young_state <- factor(result$young_state,levels = c(paste0("E",c(1:state_num))))
#   result$variable <- factor(result$variable,levels = c(paste0("E",c(1:state_num))))
#   result <- result %>%
#     arrange(young_state, variable)
#   result$label <- paste0(result$young_state,"-",result$variable)
#   result <- result[,c("label","proportion")]
#   colnames(result)[2] <- tissue_label_change(tissue)
#   if(nrow(summary)==0){
#     summary <- result
#   }else{
#     summary <- merge(summary,result,by="label")
#   }
# }
# 
# combinations <- expand.grid(paste0("E",c(1:state_num)),paste0("E",c(1:state_num)))
# combinations$Var1 <- factor(combinations$Var1,paste0("E",c(1:state_num)))
# combinations$Var2 <- factor(combinations$Var2,paste0("E",c(1:state_num)))
# combinations<- combinations %>%
#   arrange(Var1, Var2)
# combinations <- apply(combinations, 1, function(x) paste(x, collapse = "-"))
# summary$label <- factor(summary$label,levels = combinations)
# summary <- summary[order(summary$label),]
# summary <- summary[-which(summary$label %in% paste0(paste0("E",c(1:state_num)),"-",paste0("E",c(1:state_num)))),]
# rownames(summary) <- summary$label
# summary <- summary[,-1]
# summary_rowmeans <- rowMeans(summary)
# summary_rowmeans <- data.frame(condition= rownames(summary), rowmean=summary_rowmeans)
# summary_rowmeans <- summary_rowmeans[order(summary_rowmeans$rowmean),]
# summary <- summary[summary_rowmeans$condition,]
# 
# summary_colsum <- colSums(summary)
# summary_colsum <- data.frame(tissue = colnames(summary), colsum=summary_colsum)
# summary_colsum <- summary_colsum[order(summary_colsum$colsum),]
# summary <- summary[,summary_colsum$tissue]
# 
# pheatmap::pheatmap(summary,
#                    breaks = seq(0, 5, length.out = 101),
#                    fontsize = 10,na_col = "white",cluster_rows = F,cluster_cols = F,
#                    filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/all_tissue_state_transition_enrichment_score_young_replicates_heatmap_cross_comparison.png"),
#                    width = 10,height = 20)
# pheatmap::pheatmap(summary,
#                    breaks = seq(-3, 3, length.out = 101),
#                    fontsize = 10,na_col = "white",scale = "row",clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan",
#                    filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/all_tissue_state_transition_enrichment_score_young_replicates_heatmap_cross_comparison_scale_row.png"),
#                    width = 10,height = 20)
# 
# summary_part <- summary[,!(colnames(summary) %in% c("Ovary","Mammarygland","Thymus"))]
# summary_part_rowmeans <- rowMeans(summary_part)
# summary_part_rowmeans <- data.frame(condition= rownames(summary_part), rowmean=summary_part_rowmeans)
# summary_part_rowmeans <- summary_part_rowmeans[order(summary_part_rowmeans$rowmean),]
# summary_part <- summary_part[summary_part_rowmeans$condition,]
# pheatmap::pheatmap(summary_part,
#                    breaks = seq(0, 5, length.out = 101),
#                    fontsize = 10,na_col = "white",cluster_rows = F,cluster_cols = F,
#                    filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/all_tissue_state_transition_enrichment_score_young_replicates_heatmap_cross_comparison_remove_Ovary_MG_Thymus.png"),
#                    width = 10,height = 20)
# pheatmap::pheatmap(summary_part,
#                    breaks = seq(-3, 3, length.out = 101),
#                    fontsize = 10,na_col = "white",scale = "row",clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan",
#                    filename = paste0("result/all/ChromHMM/all_tissues/11_all_tissues/state_transfer/all_tissues_state_transitions/all_tissue_state_transition_enrichment_score_young_replicates_heatmap_cross_comparison_remove_Ovary_MG_Thymus_scale_row.png"),
#                    width = 10,height = 20)
# 
