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

tissues <- sort(c("mammarygland","BAT","CB","lung","kidney","aorta","brain",
                  "spleen","thymus","skin","bladder","bonemarrow","Hip",
                  "muscle","iWAT","jejunum","uterus","ovary","liver","ileum"))
tissue <- "kidney"
bin_size <- "10kb"
state_num <- "14"
to_plot <- data.frame()
for(tissue in tissues){
  opposite_change_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed"))
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  opposite_change_region$H3K9me3_condition <- "Out of H3K9me3 peaks"
  opposite_change_region$H3K9me3_condition[which(opposite_change_region$V4 %in% H3K9me3_peak_region$V4)] <- "In H3K9me3 peaks"
  
  H3K27me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_10kb_in_old_merge-W1000-G3000-E100.bed"))
  opposite_change_region$H3K27me3_condition <- "Out of H3K27me3 peaks"
  opposite_change_region$H3K27me3_condition[which(opposite_change_region$V4 %in% H3K27me3_peak_region$V4)] <- "In H3K27me3 peaks"
  
  opposite_change_region$condition <- "Neither in H3K9me3 nor in H3K27me3 peaks"
  opposite_change_region$condition[which(opposite_change_region$H3K9me3_condition=="In H3K9me3 peaks" & opposite_change_region$H3K27me3_condition=="In H3K27me3 peaks")] <- "Both in H3K9me3 and H3K27me3 peaks"
  opposite_change_region$condition[which(opposite_change_region$H3K9me3_condition=="In H3K9me3 peaks" & opposite_change_region$H3K27me3_condition=="Out of H3K27me3 peaks")] <- "In H3K9me3 peaks"
  opposite_change_region$condition[which(opposite_change_region$H3K9me3_condition=="Out of H3K9me3 peaks" & opposite_change_region$H3K27me3_condition=="In H3K27me3 peaks")] <- "In H3K27me3 peaks"
  
  t_to_plot <- as.data.frame(table(opposite_change_region$condition))
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq)*100
  t_to_plot$tissue <- tissue_label_change(tissue)
  to_plot <- rbind(to_plot,t_to_plot)
}
to_plot$Var1 <- factor(to_plot$Var1,levels=c("Both in H3K9me3 and H3K27me3 peaks","In H3K9me3 peaks","In H3K27me3 peaks", "Neither in H3K9me3 nor in H3K27me3 peaks"))
to_plot_out <- to_plot[which(to_plot$Var1=="Neither in H3K9me3 nor in H3K27me3 peaks"),]
to_plot_out <- to_plot_out[order(to_plot_out$percent),]
to_plot$tissue <- factor(to_plot$tissue,levels=to_plot_out$tissue)
color <- setNames(c("#009980","#739940","#E69900","#838B8B"),c("Both in H3K9me3 and H3K27me3 peaks","In H3K9me3 peaks","In H3K27me3 peaks", "Neither in H3K9me3 nor in H3K27me3 peaks"))
freq_sum <- to_plot %>%  
  group_by(tissue) %>%  
  summarise(Freq_sum = sum(Freq))  
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")+
  geom_text(data = freq_sum, aes(x = tissue, y = -5, label = Freq_sum),   
            inherit.aes = FALSE,  # 避免继承主 ggplot 中的 aes 映射  
            vjust = 1, size = 5, color = "black")  +
  ggtitle(paste0("H3K27me3 increase and H3K9me3 decrease bins"))

tissues <- sort(c("mammarygland","BAT","CB","lung","kidney","aorta","brain",
                  "spleen","thymus","skin","bladder","bonemarrow","Hip",
                  "muscle","iWAT","jejunum","uterus","ovary","liver"))
for(tissue in tissues){
  H3K27me3<-read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3<-read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  logFC<-merge(H3K27me3[,c("Geneid","LogFC.old.young","FDR.old.young")],H3K9me3[,c("Geneid","LogFC.old.young","FDR.old.young")],by="Geneid")
  colnames(logFC)[2:5]<-c("logFC_H3K27me3","FDR_H3K27me3","logFC_H3K9me3","FDR_H3K9me3")
  
  logFC<-merge(logFC,H3K27me3[,c("Geneid","Chr","Start","End")],by="Geneid")
  logFC<-logFC[which(logFC$FDR_H3K27me3<0.05 & logFC$FDR_H3K9me3<0.05),]
  logFC<-logFC[which(logFC$logFC_H3K27me3>0 & logFC$logFC_H3K9me3<0),]
  
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  H3K27me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  logFC<-logFC[-which(logFC$Geneid %in% H3K27me3_peak_region$V4),]
  logFC<-logFC[-which(logFC$Geneid %in% H3K9me3_peak_region$V4),]
  if(nrow(logFC) > 100){
    logFC_region <- logFC[,c("Chr","Start","End")]
    logFC_region$Start <- logFC_region$Start+1
    
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
    
    logFC_region <- as.data.table(logFC_region)
    setDT(logFC_region)
    setkey(logFC_region,Chr,Start,End)
    overlaps_young <- foverlaps(logFC_region, chromHMM_young, type = "any", nomatch = 0L)  
    overlaps_old <- foverlaps(logFC_region, chromHMM_old, type = "any", nomatch = 0L)  
    
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
      ggtitle(paste0(tissue_label_change(tissue)," H3K27me3 increase H3K9me3 decrease"), "Neither in H3K27me3 peaks nor in H3K9me3 peaks")
    print(p)
  }
}

