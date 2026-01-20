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
antibodys <- c("H3K9me3")
window_size="5000"
gap_size="10000"

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion_diff_after_remove_batch_effect.csv"))
      diff <- diff[which(diff$Length > 200000),]
    }else{
      diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_young_old_narrowpeak_diff_after_remove_batch_effect.csv"))
    }
    sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
    t_sig<-as.data.frame(table(diff$Significant))
    sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
    sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
    colnames(sig)[2] <- "Freq"
    sig <- sig[, -3] 
    
    sig$tissue <- tissue_label_change(tissue)
    sig$antibody <- antibody
    diff_peak_number<-rbind(diff_peak_number,sig)
  }
}

diff_peak_number$Freq[which(diff_peak_number$Var1=="Down")] <- -diff_peak_number$Freq[which(diff_peak_number$Var1=="Down")]
diff_peak_number <- diff_peak_number[which(diff_peak_number$Var1 %in% c("Up","Down")),]
diff_peak_number_rank <- diff_peak_number[which(diff_peak_number$Var1=="Down"),]
diff_peak_number_rank <- diff_peak_number_rank[order(diff_peak_number_rank$Freq, decreasing = T),]
diff_peak_number$tissue <- factor(diff_peak_number$tissue, levels=diff_peak_number_rank$tissue)
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(diff_peak_number_rank$tissue)))
color <- setNames(c("#e64b35","#3c5488"),c("Up","Down"))
p <- ggplot(diff_peak_number, aes(x = Freq, y = tissue, fill = Var1)) +  
  geom_bar(stat = "identity") +  
  labs(x = "Count" , y = "Tissue") +  
  theme_bw() +
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
  # geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-350, 350),breaks = seq(-350, 350, by = 100), labels = function(x) abs(x))   
p
ggsave("result/figures/H3K9me3_diff_peaks_count.pdf",p,width = 10,height = 10)  
