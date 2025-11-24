rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(ggrepel)
library(ggalluvial)  
library(scales)  
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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue <- "kidney"
state_num <- "14"
regions <- list()
regions[["second"]] <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed"))
regions[["third"]] <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_third_quadrant_after_remove_batch_effect.bed"))

H3K9me3_H3K27me3_four_quadrant_chromHMM_annotation <- function(tissue,condition,state_num){
  region <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_",condition,"_quadrant_after_remove_batch_effect.bed"))
  if(nrow(region) > 100){
    region$V2 <- region$V2+1
    region <- as.data.table(region)
    setDT(region)
    setkey(region,V1,V2,V3)
    
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
    
    overlaps_young <- foverlaps(region, chromHMM_young, type = "any", nomatch = 0L)  
    overlaps_old <- foverlaps(region, chromHMM_old, type = "any", nomatch = 0L)  
    
    overlaps_young$label <- paste0(overlaps_young$V1,":",overlaps_young$V2,"-",overlaps_young$V3)
    overlaps_old$label <- paste0(overlaps_old$V1,":",overlaps_old$V2,"-",overlaps_old$V3)
    overlaps <- merge(overlaps_young[,c("V4","label")],overlaps_old[,c("V4","label")],by="label")
    overlaps <- as.data.frame(overlaps)
    colnames(overlaps)[2:3] <- c("Young_state","Old_state") 
    
    to_plot <- overlaps %>%  
      group_by(Young_state, Old_state) %>%  
      summarize(freq = n())  
    to_plot$Young_state <- factor(to_plot$Young_state,levels=paste0("E",1:state_num))
    to_plot$Old_state <- factor(to_plot$Old_state,levels=paste0("E",1:state_num))
    colors <- hue_pal()(14) 
    colors <- setNames(colors,paste0("E",c(1:state_num)))
    
    p <- ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
      geom_alluvium(aes(fill = Young_state)) +  
      scale_fill_manual(values = colors) +
      geom_stratum() +  
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
      theme_minimal() +  
      labs(y = "Count", x = "State Transition", 
           fill = "Young State")+
      ggtitle(paste0(tissue_label_change(tissue)," ",condition," quadrant"))
    return(p)
  }
}
condition <- "second"
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- H3K9me3_H3K27me3_four_quadrant_chromHMM_annotation(tissue,condition,state_num)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list,no_of_cols = 7,no_of_rows = ceiling(length(p_list)/7))
ggsave(paste0("result/all/H3K27me3_H3K9me3/all_tissues_",condition,"_quadrant_after_remove_batch_effect_chromHMM_annotation.png"),width = 35,height = ceiling(length(p_list)/7)*5, type="cairo")



