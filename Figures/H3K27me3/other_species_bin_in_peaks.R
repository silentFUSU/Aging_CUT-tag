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
diff_peak_number <- data.frame()
for(tissue in c("MEF","BJ","MEF_EZH2_inhibit")){
  bin_size <- "10kb"
  if(tissue=="MEF"){
    peaks <- read.table(paste0("data/samples/MEF/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
    diff <- read.csv(paste0("data/samples/MEF/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  }else if(tissue=="worm"){
    peaks <- read.table(paste0("~/projects/worm/data/worm/H3K27me3/bed/H3K27me3_young_merge-W1000-G3000-E100.bed"))
    diff <- read.csv(paste0("~/projects/worm/data/worm/H3K27me3/H3K27me3_1kb_bins_diff_after_remove_batch_effect.csv"))
  }else if(tissue=="drosophila"){
    peaks <- read.table(paste0("~/projects/drosophila/data/H3K27me3/bed/H3K27me3_young_merge-W1000-G3000-E100.bed"))
    diff <- read.csv(paste0("~/projects/drosophila/data/H3K27me3/H3K27me3_1kb_bins_diff_after_remove_batch_effect.csv"))
  }else if(tissue=="BJ"){
    peaks <- read.table(paste0("~/projects/BJ_cell/data/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
    diff <- read.csv(paste0("~/projects/BJ_cell/data/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  }else if(tissue=="drosophila_female"){
    peaks <- read.table(paste0("~/projects/drosophila/data/H3K27me3_f/bed/H3K27me3_f_young_merge-W1000-G3000-E100.bed"))
    diff <- read.csv(paste0("~/projects/drosophila/data/H3K27me3_f/H3K27me3_f_1kb_bins_diff_after_remove_batch_effect.csv"))
  }else if(tissue=="BJ_bleomycin"){
    peaks <- read.table(paste0("~/projects/Aging_CUT_Tag/data/public_data/cellular_aging_GSE133292/Chip_seq/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
    diff <- read.csv("data/public_data/cellular_aging_GSE133292/Chip_seq/H3K27me3_10kb_bins_diff.csv")
  }else if(tissue=="MEF_EZH2_inhibit"){
    peaks <- read.table(paste0("~/projects/Aging_CUT_Tag/data/samples/MEF_EZH2_inhibit/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
    diff <- read.csv("~/projects/Aging_CUT_Tag/data/samples/MEF_EZH2_inhibit/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv")
  }
  
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  
  bin <- as.data.table(diff[,c("Chr","Start","End")])
  setDT(bin)
  setkey(bin,Chr,Start,End)
  
  
  overlaps <- as.data.frame(foverlaps(peaks,bin, type = "any", nomatch = 0L))
  diff_in_peaks <- diff[which(paste(diff$Chr,diff$Start,diff$End,sep = "-") %in% paste(overlaps$V1,overlaps$Start,overlaps$End,sep = "-")),]
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
  sig$antibody <- "H3K27me3"
  sig <- merge(sig,sig_in_peaks,by="Var1")
  diff_peak_number<-rbind(diff_peak_number,sig)
}

diff_peak_number$Freq[which(diff_peak_number$Var1=="Down")] <- -diff_peak_number$Freq[which(diff_peak_number$Var1=="Down")]
diff_peak_number$Freq_in_peaks[which(diff_peak_number$Var1=="Down")] <- -diff_peak_number$Freq_in_peaks[which(diff_peak_number$Var1=="Down")]
diff_peak_number <- diff_peak_number[which(diff_peak_number$Var1 %in% c("Up","Down")),]

diff_peak_number_rank <- diff_peak_number[which(diff_peak_number$Var1=="Down"),]
diff_peak_number_rank <- diff_peak_number_rank[order(diff_peak_number_rank$Freq_in_peaks),]
diff_peak_number$tissue <- factor(diff_peak_number$tissue, levels=rev(diff_peak_number_rank$tissue))
color <- setNames(c("#e64b35","#3c5488"),c("Up","Down"))
p <- ggplot(diff_peak_number, aes(x = Freq, y = tissue, fill = Var1)) +
  geom_bar(stat = "identity") +  
  geom_bar(aes(x = Freq_in_peaks), stat = "identity", fill = "black", alpha = 0.5) + 
  labs(x = "Count" , y = "Tissue") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  ) + 
  scale_fill_manual(values = color, breaks = sort(diff_peak_number_rank$tissue)) +
  geom_vline(xintercept = 0, color = "white") 
  # scale_x_continuous(limits = c(-60000, 60000),
  #                    breaks = seq(-60000, 60000, by = 30000), 
  #                    labels = function(x) format(abs(x), scientific = FALSE))   
ggsave("result/figures/H3K27me3_diff_bins_count_cellular_senescence.pdf",p,height = 3,width = 6)
