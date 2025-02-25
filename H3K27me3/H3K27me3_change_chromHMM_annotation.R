rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggalluvial)  
library(data.table)
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
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "CB"
state_num <- "14"
H3K27me3_change_region_chromHMM_annotation <- function(tissue,state_num){
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
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
  
  conditions <- c("Up","Down")
  p_list <- list()
  for(condition in conditions){
    H3K27me3_condition <- H3K27me3[which(H3K27me3$Significant==condition),c("Chr","Start","End")]
    H3K27me3_condition$Start <- H3K27me3_condition$Start+1
    H3K27me3_condition <- as.data.table(H3K27me3_condition)
    setDT(H3K27me3_condition)  
    setkey(H3K27me3_condition, Chr, Start, End) 
    
    overlaps_young <- foverlaps(H3K27me3_condition, chromHMM_young, type = "any", nomatch = 0L)  
    overlaps_old <- foverlaps(H3K27me3_condition, chromHMM_old, type = "any", nomatch = 0L)  
    
    overlaps_young$label <- paste(overlaps_young$Chr,overlaps_young$V2,overlaps_young$V3,sep = "-")
    overlaps_old$label <- paste(overlaps_old$Chr,overlaps_old$V2,overlaps_old$V3,sep = "-")
    overlaps <- merge(overlaps_young[,c("V4","label")],overlaps_old[,c("V4","label")],by="label")
    overlaps <- as.data.frame(overlaps)
    colnames(overlaps)[2:3] <- c("Young_state","Old_state") 
    to_plot <- overlaps %>%  
      group_by(Young_state, Old_state) %>%  
      summarize(freq = n())  
    to_plot$Young_state <- factor(to_plot$Young_state,levels=paste0("E",1:state_num))
    to_plot$Old_state <- factor(to_plot$Old_state,levels=paste0("E",1:state_num))
    
    p_list[[condition]] <- ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
      geom_alluvium(aes(fill = Young_state)) +  
      geom_stratum() +  
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
      theme_minimal() +  
      labs(y = "Count", x = "State Transition", 
           fill = "Young State")+
      ggtitle(paste0(tissue_label_change(tissue)," ",condition," regions"),"Sankey Plot of State Transitions")
  }
  return(p_list)
}

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
increase_p_list <- list()
decrease_p_list <- list()
for(tissue in tissues){
  p_list <- H3K27me3_change_region_chromHMM_annotation(tissue,state_num)
  increase_p_list[[tissue]] <- p_list[["Up"]]
  decrease_p_list[[tissue]] <- p_list[["Down"]]
}
increase_p_list <- increase_p_list[tissues]
decrease_p_list <- decrease_p_list[tissues]
increase_combined_plot <- plot_a_list(increase_p_list, 4, 7)
decrease_combined_plot <- plot_a_list(decrease_p_list, 4, 7)
ggsave(paste0("result/all/diff/H3K27me3/H3K27me3_increase_remove_batch_effect_chromHMM_annotation.png"),increase_combined_plot,height = 20,width = 30,type="cairo")
ggsave(paste0("result/all/diff/H3K27me3/H3K27me3_decrease_remove_batch_effect_chromHMM_annotation.png"),decrease_combined_plot,height = 20,width = 30,type="cairo")
