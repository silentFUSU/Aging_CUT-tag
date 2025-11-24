rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid) 
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
diff_number <-data.frame()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  # df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/DARs_same_peak_set/",tissue,"_DARs.txt"))
  # df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/DARs_tissue_specific/DAR_",tissue,".txt"))
  # df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_all_tissues_merge_diff_after_remove_batch_effect.csv"))
  colnames(df)[which(colnames(df)=="FDR.old.young")] <- "FDR"
  colnames(df)[which(colnames(df)=="LogFC.old.young")] <- "logFC"
  up <- df[which(df$logFC>0 & df$FDR< 0.05),]
  down <- df[which(df$logFC <0 & df$FDR < 0.05),]
  sig <- data.frame(Var1=c("Up","Down"),Freq=c(nrow(up),nrow(down)))
  colnames(sig)[2] <- "Freq"
  sig$tissue <- tissue
  sig$antibody <- "ATAC"
  diff_number<-rbind(diff_number,sig)
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
diff_number$tissue_label <- sapply(diff_number$tissue,tissue_label_change)
color <- setNames(color,sort(unique(diff_number$tissue_label)))
conditions <- c("Up","Down")

result <- diff_number %>%
  group_by(tissue_label) %>%
  summarise(Freq_sum = sum(Freq, na.rm = TRUE))
result <- result[order(result$Freq_sum),]
p_list <- list()
for(i in c(1:length(conditions))){
  condition <- conditions[i]
  df <- diff_number[which(diff_number$Var1==condition),]
  df$tissue_label <- factor(df$tissue_label,levels=result$tissue_label)
  if(i == 1){
    p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+ylab("")+
      xlab(condition)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none"  
      ) +
      scale_fill_manual(values = color) +
      theme(legend.position = "none") + 
      geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) +
      xlim(0,25000)
  }else{
    p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
      geom_bar(stat = "identity", position = position_dodge2())+
      theme_bw()+ylab("")+
      xlab(condition)+
      theme(  
        text = element_text(size = 10),  
        axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.position = "none"  
      ) + 
      scale_fill_manual(values = color) +    
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none") +   
      geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) +
      xlim(0,25000)
  }
  combined_plot <- arrangeGrob(  
    grobs = p_list,  
    ncol = length(p_list),  
    widths = c(1.4,rep(1,(length(p_list)-1))),
    top = textGrob(paste0("Differentially accessible regions number"), gp = gpar(fontsize = 15, fontface = "bold"))  
  )  
  grid.draw(combined_plot) 
}
