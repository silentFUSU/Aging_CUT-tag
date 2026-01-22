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
antibodys <- c("H3K27me3")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      bin_size <- "10kb"
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_merge-W5000-G10000-E100.bed"))
    }else{
      bin_size <- "1kb"
      peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_young_narrowpeak.bed"))
    }
    peaks <- as.data.table(peaks)
    setDT(peaks)
    setkey(peaks,V1,V2,V3)
    
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
    bin <- as.data.table(diff[,c("Chr","Start","End")])
    setDT(bin)
    setkey(bin,Chr,Start,End)
    
    overlaps <- as.data.frame(foverlaps(peaks,bin, type = "any", nomatch = 0L))
    
    # diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_fixed_bcv.csv"))
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
    
    sig$tissue <- tissue_label_change(tissue)
    sig$antibody <- antibody
    sig <- merge(sig,sig_in_peaks,by="Var1")
    diff_peak_number<-rbind(diff_peak_number,sig)
  }
}

diff_peak_number$Freq[which(diff_peak_number$Var1=="Down")] <- -diff_peak_number$Freq[which(diff_peak_number$Var1=="Down")]
diff_peak_number$Freq_in_peaks[which(diff_peak_number$Var1=="Down")] <- -diff_peak_number$Freq_in_peaks[which(diff_peak_number$Var1=="Down")]
diff_peak_number <- diff_peak_number[which(diff_peak_number$Var1 %in% c("Up","Down")),]
# 
# diff_peak_number_rank <- diff_peak_number %>%
#   group_by(tissue) %>%
#   summarize(sum_abs_freq = sum(abs(Freq)))
# 
# diff_peak_number_rank <- diff_peak_number_rank[order(diff_peak_number_rank$sum_abs_freq, decreasing = T),]

diff_peak_number_rank <- diff_peak_number[which(diff_peak_number$Var1=="Down"),]
diff_peak_number_rank <- diff_peak_number_rank[order(diff_peak_number_rank$Freq_in_peaks),]
diff_peak_number$tissue <- factor(diff_peak_number$tissue, levels=rev(diff_peak_number_rank$tissue))
color <- setNames(c("#e64b35","#3c5488"),c("Up","Down"))

p <- ggplot(diff_peak_number, aes(x = Freq, y = tissue, fill = Var1)) +
  geom_bar(stat = "identity") +  
  # geom_bar(aes(x = Freq_in_peaks), stat = "identity", fill = "black", alpha = 0.5) + 
  geom_bar_pattern(aes(x = Freq_in_peaks), 
                   stat = "identity", 
                   pattern = "stripe",              # 使用条纹图案
                   pattern_density = 0.05,          # 调整图案的密度
                   pattern_fill = NA,               # 图案内部不填充颜色
                   pattern_color = "black",         # 网格线条的颜色
                   pattern_angle = 45,              # 线条的角度
                   fill = "black", 
                   alpha = 0.5) + 
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
  geom_vline(xintercept = 0, color = "white") +
  scale_x_continuous(limits = c(-100000, 100000),
                     breaks = seq(-100000, 100000, by = 25000), 
                     labels = function(x) format(abs(x), scientific = FALSE))   
ggsave("result/figures/H3K27me3_diff_bins_count_in_young_peak.pdf",p,width = 10,height = 10)  
