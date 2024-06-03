rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(maditr)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

tissue <- "brain"
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin")
young <- c("young1","young2")
old <- c("old1","old2")
# transfer_matrix <- data.frame(tissue = character(),  
#                               young_sample = character(),
#                               young_state = character(),  
#                               old_sample = character(),
#                               old_state = character(),  
#                               Freq = numeric(),  
#                               stringsAsFactors = FALSE)  
custom_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#C0C0C0", "#808080", "#FFA500", "#FF1493")
plist <- list()
count=1
dir.create("result/all/ChromHMM/15_until_skin/state_transfer/")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  transfer_matrix <- data.frame(tissue = character(),  
                                young_sample = character(),
                                young_state = character(),  
                                old_sample = character(),
                                old_state = character(),  
                                Freq = numeric(),  
                                stringsAsFactors = FALSE)  
  for(j in c(1:length(young))){
    young_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_",young[j],"_15_segments_1k.bed"),header = F)
    states <- as.data.frame(table(young_bed$V4)) 
    states <- as.vector(states$Var1)
    for(k in c(1:length(old))){
      old_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_",old[k],"_15_segments_1k.bed"),header = F)
      for(state in states){
        state_region <- young_bed[which(young_bed$V4 == state),]
        transfer_state <- old_bed[which(paste0(old_bed$V1,"_",old_bed$V2,"_",old_bed$V3) %in% paste0(state_region$V1,"_",state_region$V2,"_",state_region$V3)),]
        t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
        colnames(t_transfer_matrix) <- c("old_state","Freq")
        t_transfer_matrix <- data.frame(tissue = tissue,young_sample = young[j], 
                                        young_state = state, old_sample = old[k], 
                                        t_transfer_matrix)
        transfer_matrix <- rbind(transfer_matrix,t_transfer_matrix)
      }
      transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",1:15)))
      transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",1:15)))
    }
  }
  result_transfer_matrix <- transfer_matrix[which(transfer_matrix$young_state!=transfer_matrix$old_state),] %>%  
    group_by(young_state, old_state) %>%  
    summarise(mean_freq = mean(Freq, na.rm = TRUE))  
  result_transfer_matrix2 <- result_transfer_matrix %>%  
    group_by(young_state) %>%  
    mutate(percent_freq = mean_freq/sum(mean_freq)) 
  plist[[count]] <- ggplot(result_transfer_matrix2, aes(x = young_state, weight = percent_freq, fill = old_state))+
    geom_bar( position = "stack")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values =custom_colors)+ggtitle(tissue)
  count=count+1
}

dcombined_plot <- plot_a_list(plist, 4, 5)
ggsave(paste0("result/all/ChromHMM/15_until_skin/state_transfer/all_tissue_state_transfer_only_transfer.png"),width = 40,height = 40,type="cairo")

transfer_matrix <- data.frame(tissue = character(),  
                              young_sample = character(),
                              young_state = character(),  
                              old_sample = character(),
                              old_state = character(),  
                              Freq = numeric(),  
                              stringsAsFactors = FALSE)  
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(young))){
    young_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_",young[j],"_15_segments_1k.bed"),header = F)
    states <- as.data.frame(table(young_bed$V4)) 
    states <- as.vector(states$Var1)
    for(k in c(1:length(old))){
      old_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_",old[k],"_15_segments_1k.bed"),header = F)
      for(state in states){
        state_region <- young_bed[which(young_bed$V4 == state),]
        transfer_state <- old_bed[which(paste0(old_bed$V1,"_",old_bed$V2,"_",old_bed$V3) %in% paste0(state_region$V1,"_",state_region$V2,"_",state_region$V3)),]
        t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
        colnames(t_transfer_matrix) <- c("old_state","Freq")
        t_transfer_matrix <- data.frame(tissue = tissue,young_sample = young[j], 
                                        young_state = state, old_sample = old[k], 
                                        t_transfer_matrix)
        transfer_matrix <- rbind(transfer_matrix,t_transfer_matrix)
      }
      transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",1:15)))
      transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",1:15)))
    }
  }
}
# write.csv(transfer_matrix,"result/all/ChromHMM/15_until_skin/state_transfer/state_transfer.csv",row.names = F)
result_transfer_matrix <- transfer_matrix %>%
  group_by(young_state, old_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_transfer_matrix2 <- result_transfer_matrix %>%
  group_by(young_state) %>%
  mutate(percent_freq = mean_freq/sum(mean_freq))
result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
labels <-as.character(result_transfer_matrix_to_plot$young_state)
result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
pheatmap::pheatmap(result_transfer_matrix_to_plot,cluster_rows = F,cluster_cols = F, breaks = seq(0, 0.2, length.out = 101),display_numbers = T,labels_row = labels,labels_col = labels,fontsize = 10)

transfer_matrix <- read.csv("result/all/ChromHMM/15_until_skin/state_transfer/state_transfer.csv")
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",1:15)))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",1:15)))
dir.create(paste0("result/all/ChromHMM/15_until_skin/state_transfer/each_tissue_transfer/"))
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue),]
  result_transfer_matrix <- t_transfer_matrix %>%
    group_by(young_state, old_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_transfer_matrix2 <- result_transfer_matrix %>%
    group_by(young_state) %>%
    mutate(percent_freq = mean_freq/sum(mean_freq))
  result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
  labels <-as.character(result_transfer_matrix_to_plot$young_state)
  result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
  if(tissue=="brain") tissue <- "FC"
  p<-pheatmap::pheatmap(result_transfer_matrix_to_plot,cluster_rows = F,
                     cluster_cols = F, breaks = seq(0, 0.2, length.out = 101),
                     display_numbers = T,labels_row = labels,labels_col = labels,
                     fontsize = 10,main=tissue,fontsize_number = 10)
  save_pheatmap_pdf(p,paste0("result/all/ChromHMM/15_until_skin/state_transfer/each_tissue_transfer/",tissue,"_state_transfer.pdf"),width=7, height=7)
  }
