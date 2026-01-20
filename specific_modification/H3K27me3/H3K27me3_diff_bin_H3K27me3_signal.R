rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
library(dplyr)
library(GenomeInfoDb)
library("GenomicRanges")
library(genomation)
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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
for(tissue in tissues){
  diff <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  
  bin <- as.data.table(read.table("~/ref_data/mm10_10kb_bins.bed"))
  setDT(bin)
  setkey(bin,V1,V2,V3)
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peak <- peak[,c(1:3)]
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  overlaps <- foverlaps(peak, bin, type = "any", nomatch = 0L)  
  
  diff <- diff[which(diff$Geneid %in% overlaps$V4),c("Geneid","logCPM","LogFC.old.young","Significant")]
  diff$tissue <- tissue_label_change(tissue)  
  tissue_summary <- rbind(tissue_summary,diff)
    
  }
to_plot <- tissue_summary
to_plot$Significant <- factor(to_plot$Significant,levels=c("Up","Stable","Down"))
p <- ggplot(to_plot, aes(x = Significant, y = logCPM,fill=Significant)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "CPM", title = "H3K27me3 signal in diff bins") +
  theme_bw()
p

p <- ggplot(to_plot,aes(x=logCPM,y=LogFC.old.young))+    
  geom_jitter(size = 3, alpha = 0.7)+
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("logFC")+
  ylab("")+
  labs(fill = "", color = "")


p <- ggplot(to_plot[which(to_plot$Significant !="Stable"),],aes(x=logCPM,y=LogFC.old.young,color=Significant))+    
  geom_jitter(size = 3, alpha = 0.7)+
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("logCPM")+
  ylab("logFC")+
  labs(fill = "", color = "")
p

smoothScatter(to_plot$LogFC.old.young ~ to_plot$logCPM,
              nrpoints = 1000)

to_plot <- tissue_summary %>%
  group_by(tissue, Significant) %>%
  summarise(
    avg_logCPM = mean(logCPM, na.rm = TRUE),
    avg_LogFC = mean(LogFC.old.young, na.rm = TRUE)
  )

to_plot$Significant <- factor(to_plot$Significant,levels=c("Up","Stable","Down"))
p <- ggplot(to_plot,aes(x=avg_logCPM,y=avg_LogFC,color=Significant))+    
  geom_jitter(size = 3, alpha = 0.7)+
  geom_smooth(data = to_plot, aes(x = avg_logCPM, y = avg_LogFC),
              method = "lm", color = "#e64b35", se = TRUE, level = 0.95) +
  theme_bw()+theme(text = element_text(size = 18))+
  xlab("logCPM")+
  ylab("logFC")+
  labs(fill = "", color = "")
p

### only consider young signal
