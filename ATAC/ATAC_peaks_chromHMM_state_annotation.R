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
library(edgeR)
library(ggalluvial)  
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
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
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "lung"
state_num <- 14
ATAC_chromHMM_state_annotation <- function(tissue,state_num){
  peaks <- read.table("data/samples/ATAC/ATAC_peak_from_MJ/All_Samples.fwp.filter.non_overlapping.bed")
  # peaks <- peaks[grepl(tissue, peaks$V4), ]
  peaks <- peaks[,c(1:3)]
  peaks$V2 <- peaks$V2 + 1
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues_previous/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues_previous/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  
  setDT(chromHMM_young) 
  setkey(chromHMM_young, V1, V2, V3) 
  setDT(chromHMM_old) 
  setkey(chromHMM_old, V1, V2, V3) 
  
  overlaps_young <- foverlaps(peaks, chromHMM_young, type = "any", nomatch = 0L)  
  overlaps_old <- foverlaps(peaks, chromHMM_old, type = "any", nomatch = 0L)  
  
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
  
  p <- ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
    geom_alluvium(aes(fill = Young_state)) +  
    geom_stratum() +  
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
    theme_minimal() +  
    labs(y = "Count", x = "State Transition", 
         fill = "Young State")+
    ggtitle(paste0(tissue_label_change(tissue)," ATAC peaks"),"Sankey Plot of State Transitions")
  return(p)
}

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
state_num <- 14
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- ATAC_chromHMM_state_annotation(tissue,state_num)    
}
combined_plot <- plot_a_list(p_list, 4, 7)
ggsave(paste0("result/all/diff/ATAC/ATAC_Union_peaks_chromHMM_annotation.png"),combined_plot,height = 20,width = 25,type="cairo")


ATAC_change_chromHMM_annotation <- function(tissue,condition,state){
  peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_MJ/Peak_",condition,"_",tissue,".bed"))
  # peaks <- peaks[grepl(tissue, peaks$V4), ]
  # peaks <- peaks[,c(1:3)]
  peaks$V2 <- peaks$V2 + 1
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues_previous/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues_previous/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  
  setDT(chromHMM_young) 
  setkey(chromHMM_young, V1, V2, V3) 
  setDT(chromHMM_old) 
  setkey(chromHMM_old, V1, V2, V3) 
  
  overlaps_young <- foverlaps(peaks, chromHMM_young, type = "any", nomatch = 0L)  
  overlaps_old <- foverlaps(peaks, chromHMM_old, type = "any", nomatch = 0L)  
  
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
  
  p <- ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
    geom_alluvium(aes(fill = Young_state)) +  
    geom_stratum() +  
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
    theme_minimal() +  
    labs(y = "Count", x = "State Transition", 
         fill = "Young State")+
    ggtitle(paste0(tissue_label_change(tissue)," ",condition," ATAC peaks"),"Sankey Plot of State Transitions")
  return(p)
}



conditions <- c("Up","Down")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
state_num <- 14
p_list <- list()
for(condition in conditions){
  for(tissue in tissues){
    p_list[[tissue]] <- ATAC_change_chromHMM_annotation(tissue,condition,state_num)    
  }
  combined_plot <- plot_a_list(p_list, 4, 7)
  if(condition == "Up"){
    condition_label <- "increase"
  }else{
    condition_label <- "decrease"
  }
  ggsave(paste0("result/all/diff/ATAC/ATAC_",condition_label,"_peaks_remove_batch_effect_chromHMM_annotation.png"),combined_plot,height = 20,width = 25,type="cairo")
  
}