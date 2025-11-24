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
diff_peak_number <-data.frame()
mm10_10k <- read.delim("~/ref_data/mm10_10kb_bins.bed")
mm10_1k <- read.delim("~/ref_data/mm10_1kb_bins.bed")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      bin_size <- "10kb"
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_young_old_merge-W1000-G3000-E100.bed"))
    }else{
      bin_size <- "1kb"
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_1kb_in_young_old_merge_macs_narrowpeak.bed"))
    }
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
    diff_in_peaks <- diff[which(paste(diff$Chr,diff$Start,diff$End,sep = "-") %in% paste(peaks[,1],peaks[,2],peaks[,3],sep = "-")),]
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    
    sig_in_peaks <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig_in_peaks<-as.data.frame(table(diff_in_peaks$Significant))
    sig_in_peaks <- merge(sig_in_peaks, t_sig_in_peaks, by="Var1", all.x=TRUE) 
    sig_in_peaks$Freq.x <- ifelse(is.na(sig_in_peaks$Freq.y), sig_in_peaks$Freq.x, sig_in_peaks$Freq.y)  
    colnames(sig_in_peaks)[2] <- "Freq_in_peaks"
    sig_in_peaks <- sig_in_peaks[, -3] 
    
    sig$tissue <- tissue
    sig$antibody <- antibody
    sig <- merge(sig,sig_in_peaks,by="Var1")
    diff_peak_number<-rbind(diff_peak_number,sig)
  }
}
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
diff_peak_number$tissue_label <- sapply(diff_peak_number$tissue,tissue_label_change)
color <- setNames(color,sort(unique(diff_peak_number$tissue_label)))
p_list <- list()
conditions <- c("Down","Up")
diff_peak_number$peak_percent <- diff_peak_number$Freq_in_peaks/diff_peak_number$Freq*100
diff_peak_number$peak_percent_label <- paste0(diff_peak_number$Freq,"(",round(diff_peak_number$peak_percent,1),"%)")
diff_peak_number$peak_label <- paste0(diff_peak_number$Freq_in_peaks,"/",diff_peak_number$Freq)

for(condition in conditions){
  for(i in c(1:length(antibodys))){
    df <- diff_peak_number[which(diff_peak_number$antibody==antibodys[i] & diff_peak_number$Var1==condition),]
    if(i == 1){
      p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
        geom_bar(stat = "identity", position = position_dodge2())+
        theme_bw()+ylab("")+
        geom_bar(aes(x = Freq_in_peaks), stat = "identity", fill = "black", alpha = 0.5) + 
        xlab(antibodys[[i]])+
        theme(  
          text = element_text(size = 10),  
          axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
          axis.title.y = element_blank(),  
          axis.ticks.y = element_blank(),  
          legend.position = "none"  
        ) +
        scale_fill_manual(values = color) +
        theme(legend.position = "none") + 
        geom_text(aes(label = peak_label), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) +
        xlim(0,100000)
    }else{
      p_list[[i]] <-  ggplot(df,mapping = aes(x=Freq,y=tissue_label,fill = tissue_label))+
        geom_bar(stat = "identity", position = position_dodge2())+
        theme_bw()+ylab("")+
        geom_bar(aes(x = Freq_in_peaks), stat = "identity", fill = "black", alpha = 0.5) +
        xlab(antibodys[[i]])+
        theme(  
          text = element_text(size = 10),  
          axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
          axis.title.y = element_blank(),  
          axis.ticks.y = element_blank(),  
          legend.position = "none"  
        ) + scale_fill_manual(values = color) +    
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),legend.position = "none") +   
        geom_text(aes(label = peak_label), position = position_dodge2(width = 0.9), hjust = 0.1, size = 3) +
        xlim(0,100000)
    }
  }
  combined_plot <- arrangeGrob(  
    grobs = p_list,  
    ncol = length(p_list),  
    widths = c(1.4,rep(1,(length(p_list)-1))),
    top = textGrob(paste0(condition," bin number"), gp = gpar(fontsize = 15, fontface = "bold"))  
  )  
  ggsave(paste0("result/Sup_figures/all_tissues_",condition,"_after_remove_batch_effect.pdf"), plot = combined_plot, width = 18, height = 6)
}
