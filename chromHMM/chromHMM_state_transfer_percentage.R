rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(reshape2)
# library(randomcoloR)
# colors<-distinctColorPalette(20)
# write.table(colors,"data/samples/20_distinct_color.txt",row.names = F,col.names = F)
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
tissues <-  c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
              "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
plist <- list()
count=1
state_num=14
dir.create(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/"))
colors <- read.table("data/samples/20_distinct_color.txt")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  transfer_matrix <- data.frame(tissue = character(),  
                                young_state = character(),  
                                old_state = character(),  
                                Freq = numeric(),  
                                stringsAsFactors = FALSE)  
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  
  young$label <- paste0(young$V1,"-",young$V2,"-",young$V3)
  old$label <- paste0(old$V1,"-",old$V2,"-",old$V3)
  young <- young[which(young$label %in% old$label),]  
  old <- old[which(old$label %in% young$label),]
  states <- sort(unique(young$V4))
  for(state in states){
    state_region <- young[which(young$V4 == state),]
    transfer_state <- old[which(old$label %in% state_region$label),]
    t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
    colnames(t_transfer_matrix) <- c("old_state","Freq")
    t_transfer_matrix <- data.frame(tissue = tissue_label_change(tissue), 
                                    young_state = state, 
                                    t_transfer_matrix)
    transfer_matrix <- rbind(transfer_matrix,t_transfer_matrix)
  }
  transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:length(states)))))
  transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:length(states)))))
  result_transfer_matrix <- transfer_matrix[which(transfer_matrix$young_state != transfer_matrix$old_state),] %>%  
    group_by(young_state) %>%  
    mutate(percent_freq = Freq/sum(Freq)) 
  custom_colors <- colors$V1[1:length(states)]
  plist[[count]] <- ggplot(result_transfer_matrix, aes(x = young_state, weight = percent_freq, fill = old_state))+
    geom_bar( position = "stack")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 18))+
    scale_fill_manual(values =custom_colors)+ggtitle(tissue_label_change(tissue))
  count=count+1
}
combined_plot <- plot_a_list(plist, 4, 7)
ggsave(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/all_tissue_state_transfer_only_transfer.png"),combined_plot,width = 45,height = 40,type="cairo")



transfer_matrix <- data.frame(tissue = character(),  
                              young_state = character(),  
                              old_state = character(),  
                              Freq = numeric(),  
                              stringsAsFactors = FALSE)  
tissues <-   c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                         "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
state_num <- "14"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste0(young1$V1,"-",young1$V2,"-",young1$V3,"-",young1$V4)
  young2$label <- paste0(young2$V1,"-",young2$V2,"-",young2$V3,"-",young2$V4)
  young <- young1[which(young1$label %in% young2$label),]
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste0(old1$V1,"-",old1$V2,"-",old1$V3,"-",old1$V4)
  old2$label <- paste0(old2$V1,"-",old2$V2,"-",old2$V3,"-",old2$V4)
  old <- old1[which(old1$label %in% old2$label),]
  young$label <- paste0(young$V1,"-",young$V2,"-",young$V3)
  old$label <- paste0(old$V1,"-",old$V2,"-",old$V3)
  young <- young[which(young$label %in% old$label),]  
  old <- old[which(old$label %in% young$label),]
  states <- sort(unique(young$V4))
  for(state in states){
    state_region <- young[which(young$V4 == state),]
    transfer_state <- old[which(old$label %in% state_region$label),]
    t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
    colnames(t_transfer_matrix) <- c("old_state","Freq")
    t_transfer_matrix <- data.frame(tissue = tissue_label_change(tissue), 
                                    young_state = state, 
                                    t_transfer_matrix)
    transfer_matrix <- rbind(transfer_matrix,t_transfer_matrix)
  }
  transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",c(1:length(states)))))
  transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",c(1:length(states)))))
}

# write.csv(transfer_matrix,paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"),row.names = F)

#scale by row
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
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

result_transfer_matrix <- transfer_matrix[which(transfer_matrix$young_state != transfer_matrix$old_state),] %>%
  group_by(young_state, old_state) %>%
  summarise(mean_freq = mean(Freq, na.rm = TRUE))
result_transfer_matrix2 <- result_transfer_matrix %>%
  group_by(young_state) %>%
  mutate(percent_freq = mean_freq/sum(mean_freq))
result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq") 
labels <-as.character(result_transfer_matrix_to_plot$young_state)
result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
result_transfer_matrix_to_plot[is.na(result_transfer_matrix_to_plot)] <- 0
for(i in 1:nrow(result_transfer_matrix_to_plot)){
  result_transfer_matrix_to_plot[i,i] <- NA
}
pheatmap::pheatmap(result_transfer_matrix_to_plot,
                   cluster_rows = F,cluster_cols = F, 
                   breaks = seq(0, 1, length.out = 101),
                   display_numbers = T,labels_row = labels,
                   labels_col = labels,fontsize = 10,na_col = "white")

transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",1:state_num)))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",1:state_num)))
dir.create(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/each_tissue_transfer/"))
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue_label_change(tissue)),]
  result_transfer_matrix <- t_transfer_matrix %>%
    group_by(young_state, old_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_transfer_matrix2 <- result_transfer_matrix %>%
    group_by(young_state) %>%
    mutate(percent_freq = mean_freq/sum(mean_freq))
  result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq")
  labels <-as.character(result_transfer_matrix_to_plot$young_state)
  result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
  result_transfer_matrix_to_plot[is.na(result_transfer_matrix_to_plot)] <- 0
  p<-pheatmap::pheatmap(result_transfer_matrix_to_plot,cluster_rows = F,
                     cluster_cols = F, breaks = seq(0, 0.2, length.out = 101),
                     display_numbers = T,labels_row = labels,labels_col = labels,
                     fontsize = 10,main=tissue_label_change(tissue),fontsize_number = 10)
  save_pheatmap_pdf(p,paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/each_tissue_transfer/",tissue,"_state_transfer.pdf"),width=7, height=7)
  }

#scale by all transfer condition
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
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

transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",1:state_num)))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",1:state_num)))
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue_label_change(tissue)),]
  result_transfer_matrix <- t_transfer_matrix %>%
    group_by(young_state, old_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_transfer_matrix2 <- result_transfer_matrix 
  result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
  result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq")
  labels <-as.character(result_transfer_matrix_to_plot$young_state)
  result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
  result_transfer_matrix_to_plot[is.na(result_transfer_matrix_to_plot)] <- 0
  p<-pheatmap::pheatmap(result_transfer_matrix_to_plot,cluster_rows = F,
                        cluster_cols = F, breaks = seq(0, 5, length.out = 101),
                        display_numbers = T,labels_row = labels,labels_col = labels,
                        fontsize = 10,main=tissue_label_change(tissue),fontsize_number = 10)
  dir.create(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/each_tissue_transfer_scale_all_condition/"),showWarnings = F,recursive = T)
  save_pheatmap_pdf(p,paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/each_tissue_transfer_scale_all_condition/",tissue,"_state_transfer.pdf"),width=7, height=7)
}

### remove empty state
transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix <- transfer_matrix[-which(transfer_matrix$young_state %in% c("E3","E12") | transfer_matrix$old_state %in% c("E3","E12")),]
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
to_plot <- result_transfer_matrix2[which(result_transfer_matrix2$young_state != result_transfer_matrix2$old_state),]
to_plot$label <- paste(to_plot$young_state,to_plot$old_state,sep = "-")
to_plot$label <- factor(to_plot$label,levels=to_plot$label)
ggplot(to_plot, aes(x = label, y = percent_freq,color=label)) +  
  geom_point(size = 3) +  
  labs(x = NULL, y = "Percentage(%)") +  
  theme_bw()+
  ylim(0,5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 10),legend.position = "none")


transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
transfer_matrix <- transfer_matrix[-which(transfer_matrix$young_state %in% c("E3","E12") | transfer_matrix$old_state %in% c("E3","E12")),]
transfer_matrix$old_state <- factor(transfer_matrix$old_state, levels=c(paste0("E",1:state_num)))
transfer_matrix$young_state <- factor(transfer_matrix$young_state, levels=c(paste0("E",1:state_num)))
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  t_transfer_matrix <- transfer_matrix[which(transfer_matrix$tissue==tissue_label_change(tissue)),]
  result_transfer_matrix <- t_transfer_matrix %>%
    group_by(young_state, old_state) %>%
    summarise(mean_freq = mean(Freq, na.rm = TRUE))
  result_transfer_matrix2 <- result_transfer_matrix 
  result_transfer_matrix2$percent_freq <- result_transfer_matrix2$mean_freq/sum(result_transfer_matrix2$mean_freq)*100
  result_transfer_matrix_to_plot <- dcast(result_transfer_matrix2, formula = young_state~old_state, value.var = "percent_freq")
  labels <-as.character(result_transfer_matrix_to_plot$young_state)
  result_transfer_matrix_to_plot <- result_transfer_matrix_to_plot[,-1]
  result_transfer_matrix_to_plot[is.na(result_transfer_matrix_to_plot)] <- 0
  p<-pheatmap::pheatmap(result_transfer_matrix_to_plot,cluster_rows = F,
                        cluster_cols = F, breaks = seq(0, 5, length.out = 101),
                        display_numbers = T,labels_row = labels,labels_col = labels,
                        fontsize = 10,main=tissue_label_change(tissue),fontsize_number = 10)
  dir.create(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/each_tissue_transfer_scale_all_condition/"),showWarnings = F,recursive = T)
  save_pheatmap_pdf(p,paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/each_tissue_transfer_scale_all_condition/",tissue,"_remove_empty_state_transfer.pdf"),width=7, height=7)
}
