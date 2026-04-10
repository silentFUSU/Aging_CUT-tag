rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)

search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
domain_list <- data.frame()
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K27me3"

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

peak_summary <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W5000-G10000-E100.bed"))    
  df$length <- df$V3 - df$V2 +1
  t_summary <- data.frame(tissue=tissue_label_change(tissue),count=nrow(df),coverage=sum(df$length))
  peak_summary <- rbind(peak_summary,t_summary)
}
peak_summary$label <- "Peak"

# domain_summary <- data.frame()
# for(tissue in tissues){
#   domains <- read.table(paste0("data/samples/",tissue,"/",antibody,"/peaks/edd/edd_peaks_fdr05.bed"))
#   domains$length <- domains$V3 - domains$V2 +1
#   domains <- domains[which(domains$length >= 200000),]
#   t_summary <- data.frame(tissue=tissue_label_change(tissue),count=nrow(domains),coverage=sum(domains$length))
#   domain_summary <- rbind(domain_summary,t_summary)
# }
# domain_summary$label <- "Domain"
# to_plot <- rbind(peak_summary,domain_summary)
to_plot <- peak_summary
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$tissue)))
to_plot$label <- factor(to_plot$label,levels=c("Peak","Domain"))
p <- ggplot(to_plot, aes(x = label,y = count),) +
  geom_violin(fill = "lightgray", color = NA, alpha = 0.5) +
  geom_boxplot(width=0.1) +
  geom_jitter(aes(color = tissue), width = 0.2, size = 2) + 
  scale_color_manual(values = color) +
  labs(y = "count",x=NULL) +
  theme_bw()+
  theme(
  )+ylim(0,12000)
p
ggsave("result/Sup_figures/H3K27me3_peak_number.pdf",p,width = 6,height = 8)
