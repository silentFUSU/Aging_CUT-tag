rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(maditr)
state_num <- 11
tissue <- "lung"
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
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
chromHMM_state_transfer_H3K27me3_H3K9me3_change <- function(tissue,state_num){
  H3K27me3_H3K9me3 <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant.bed"))
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
  
  young$label <- paste0(young$V1,":",young$V2,"-",young$V3)
  old$label <- paste0(old$V1,":",old$V2,"-",old$V3)
  young <- young[which(young$label %in% old$label),]  
  old <- old[which(old$label %in% young$label),]
  colnames(young)[4] <- "young_state"
  colnames(old)[4] <- "old_state"
  state_change <- merge(young[,c(1:5)],old[,c(4:5)],by="label")
  state_change <- state_change[which(state_change$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
  state_change <- state_change[,c(-1)]
  state_change$condition <- paste0(state_change$young_state,"-",state_change$old_state)
  
  H3K27me3_H3K9me3$V2 <- H3K27me3_H3K9me3$V2 + 1
  H3K27me3_H3K9me3 <- as.data.table(H3K27me3_H3K9me3)
  setDT(H3K27me3_H3K9me3)  
  setkey(H3K27me3_H3K9me3,V1,V2,V3)
  
  state_change$V2 <- state_change$V2 + 1 
  state_change <- as.data.table(state_change)
  setDT(state_change) 
  setkey(state_change,V1,V2,V3)
  
  overlaps <- foverlaps(H3K27me3_H3K9me3,state_change, type = "any", nomatch = 0L)
  overlaps <- overlaps[which(overlaps$young_state != overlaps$old_state),]
  to_plot <- as.data.frame(table(overlaps$condition))
  to_plot <- to_plot[which(to_plot$Freq >100),]
  to_plot$percent <- to_plot$Freq/sum(to_plot$Freq)*100
  ggplot(to_plot, aes(x=1,y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)," (",nrow(overlaps)," regions)"))
}


sanky_plot <- function(tissue,state_num){
  histone <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed"))
  if(nrow(histone) < 200){
    return(0)
  }
  histone$V2 <- histone$V2+1
  histone <- as.data.table(histone)
  setDT(histone)  
  setkey(histone, V1, V2, V3) 
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  setDT(chromHMM_young) 
  setkey(chromHMM_young, V1, V2, V3) 
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  setDT(chromHMM_old) 
  setkey(chromHMM_old, V1, V2, V3) 
  
  overlaps_young <- foverlaps(histone, chromHMM_young, type = "any", nomatch = 0L)  
  overlaps_old <- foverlaps(histone, chromHMM_old, type = "any", nomatch = 0L)  
  
  overlaps_young$label <- paste(overlaps_young$V1,overlaps_young$V2,overlaps_young$V3,sep = "-")
  overlaps_old$label <- paste(overlaps_old$V1,overlaps_old$V2,overlaps_old$V3,sep = "-")
  overlaps <- merge(overlaps_young[,c("V4","label")],overlaps_old[,c("V4","label")],by="label")
  overlaps <- as.data.frame(overlaps)
  colnames(overlaps)[2:3] <- c("Young_state","Old_state") 
  to_plot <- overlaps %>%  
    group_by(Young_state, Old_state) %>%  
    summarize(freq = n())  
  to_plot$Young_state <- factor(to_plot$Young_state,levels=paste0("E",1:state_num))
  to_plot$Old_state <- factor(to_plot$Old_state,levels=paste0("E",1:state_num))
  color <- read.table("data/samples/20_distinct_color.txt")
  color <- setNames(color$V1,paste0("E",1:11))
  p <- ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
    geom_alluvium(aes(fill = Young_state)) +  
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
    theme_minimal() +  
    scale_fill_manual(values = color) +
    labs(y = "Count", x = "State Transition", 
         fill = "Young State")+
    theme(  
      plot.title = element_text(size = 14),          
      axis.title = element_text(size = 12),         
      axis.text = element_text(size = 10),           
      legend.title = element_text(size = 10),      
      legend.text = element_text(size = 9)           
    ) +
    ggtitle(paste0(tissue_label_change(tissue)," (",nrow(overlaps)," regions)"),"H3K27me3 increased and H3K9me3 decreased regions")
  return(p)
}

p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- sanky_plot(tissue,state_num)
}
is_not_double <- function(x) {   
  !is.double(x)  
}  
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
filtered_list <- Filter(is_not_double, p_list)  
combined_plot <- plot_a_list(filtered_list, 4, 5)
ggsave("result/all/H3K27me3_H3K9me3/all_tissue_H3K27me3_H3K9me3_after_remove_batch_effect_chromHMM_annotation.png",combined_plot,width = 25,height = 20,type="cairo")

summary <- data.frame()
for(tissue in tissues){
  H3K27me3_H3K9me3 <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed"))
  if(nrow(H3K27me3_H3K9me3) > 200){
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
    
    young$label <- paste0(young$V1,":",young$V2,"-",young$V3)
    old$label <- paste0(old$V1,":",old$V2,"-",old$V3)
    young <- young[which(young$label %in% old$label),]  
    old <- old[which(old$label %in% young$label),]
    colnames(young)[4] <- "young_state"
    colnames(old)[4] <- "old_state"
    state_change <- merge(young[,c(1:5)],old[,c(4:5)],by="label")
    state_change <- state_change[which(state_change$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    state_change <- state_change[,c(-1)]
    state_change$condition <- paste0(state_change$young_state,"-",state_change$old_state)
    
    H3K27me3_H3K9me3$V2 <- H3K27me3_H3K9me3$V2 + 1
    H3K27me3_H3K9me3 <- as.data.table(H3K27me3_H3K9me3)
    setDT(H3K27me3_H3K9me3)  
    setkey(H3K27me3_H3K9me3,V1,V2,V3)
    
    state_change$V2 <- state_change$V2 + 1 
    state_change <- as.data.table(state_change)
    setDT(state_change) 
    setkey(state_change,V1,V2,V3)
    
    overlaps <- foverlaps(H3K27me3_H3K9me3,state_change, type = "any", nomatch = 0L)
    overlaps <- overlaps[which(overlaps$young_state != overlaps$old_state),]
    if(nrow(overlaps)>0){
      to_plot <- as.data.frame(table(overlaps$condition))
      to_plot <- to_plot[which(to_plot$Freq >100),]
      if(nrow(to_plot) > 0){
        to_plot$percent <- to_plot$Freq/sum(to_plot$Freq)*100
        t_summary <- to_plot[,c("Var1","percent")]
        colnames(t_summary)[2] <- tissue_label_change(tissue)
        if(nrow(summary) ==0){
          summary <- t_summary
        }else{
          summary <- merge(summary,t_summary,by="Var1",all=T)
        }
      }
    } 
  }
}
to_plot <- reshape2::melt(summary)
to_plot_filtered <- to_plot[!is.na(to_plot[, 3]), ]  
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$Var1)))
ggplot(to_plot_filtered, aes(x=variable,y = value, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values=color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),legend.title = element_blank()) +
  ylab("Proportion")+
  ggtitle("H3K27me3 increased and H3K9me3 decreased regions","chromHMM state changed proportion")
