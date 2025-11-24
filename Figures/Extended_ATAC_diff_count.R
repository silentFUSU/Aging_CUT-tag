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
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
  up <- df[which(df$LogFC.old.young > 0 & df$FDR.old.young < 0.05),]
  down <- df[which(df$LogFC.old.young < 0 & df$FDR.old.young < 0.05),]
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
conditions <- c("Up","Down")
color <- setNames(c("#4589C8FF","#EE7C7AFF"),c("Down","Up"))
p_list <- list()
df <- diff_number[which(diff_number$Var1 %in% conditions),]
df$tissue_label <- factor(df$tissue_label,levels=rev(sort(unique(df$tissue_label))))
p <- ggplot(df,mapping =  aes(x=Freq,y=tissue_label,fill = Var1))+
  # geom_bar(stat = "identity",position = position_dodge2(),aes(alpha = ifelse(Var1 == "Down", 0.9, 1)),color = "black")+
  geom_bar(stat = "identity",position = position_dodge2())+
  theme_bw()+ylab("")+
  xlab("ATAC")+
  theme(  
    text = element_text(size = 10),  
    axis.text.y = element_text(size = 11),  # 改变 y 轴刻度标签的字体大小  
    axis.title.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    legend.position = "none"  
  ) +
  scale_fill_manual(values = color) +
  theme(legend.position = "none") +
  scale_x_continuous(
    limits = c(0,12000), 
    breaks = c(0,12000)
  )
ggsave(paste0("result/figures/all_tissues_DAR_after_remove_batch_effect.pdf"), plot = p, width = 4, height = 6)

diff_number$Freq[which(diff_number$Var1=="Down")] <- -diff_number$Freq[which(diff_number$Var1=="Down")]
diff_number_rank <- diff_number %>%
  group_by(tissue_label) %>%
  summarize(sum_abs_freq = sum(abs(Freq)))
diff_number_rank <- diff_number_rank[order(diff_number_rank$sum_abs_freq, decreasing = T),]
diff_number$tissue_label <- factor(diff_number$tissue_label, levels=rev(diff_number_rank$tissue_label))
color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(diff_number_rank$tissue_label)))

ggplot(diff_number, aes(x = Freq, y = tissue_label, fill = tissue_label)) +
  geom_bar(stat = "identity") +  
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
  scale_fill_manual(values = color)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-15000, 15000),
                     breaks = seq(-15000, 15000, by = 5000), 
                     labels = function(x) format(abs(x), scientific = FALSE)) 
                  