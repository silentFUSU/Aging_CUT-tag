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
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/10kb_bins_all_depth.csv"))
  # df <- df[which(df$total_V5 > 15),]
  bin <- as.data.table(read.table("~/ref_data/mm10_10kb_bins.bed"))
  setDT(bin)
  setkey(bin,V1,V2,V3)
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peak <- peak[,c(1:3)]
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  overlaps <- foverlaps(peak, bin, type = "any", nomatch = 0L)  
  
  df <- df[which(df$label %in% overlaps$V4),]
  
  diff <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  diff <- diff[,c("Geneid","Significant")]
  df <- merge(df,diff,by.x="label",by.y="Geneid")
  df <- merge(df,search_table,by.x="sample",by.y="sample_name")
  tissue_summary <- rbind(tissue_summary,df)
  }

to_plot <- tissue_summary[which(tissue_summary$age=="3M"),]
to_plot$Significant <- factor(to_plot$Significant,levels=c("Up","Stable","Down"))
to_plot$percent <- to_plot$percent *100
p <- ggplot(to_plot, aes(x = Significant, y = percent,fill=Significant)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "DNA methylation", title = "all tissues DNA methylation") +
  theme_bw()+
  ylim(50,100)
p
ggsave("result/Sup_figures/all_tissues_H3K27me3_young_DNA_methylation.pdf",p,width = 6,height = 8)

##delta
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/10kb_bins_all_depth.csv"))
  # df <- df[which(df$total_V5 > 15),]
  bin <- as.data.table(read.table("~/ref_data/mm10_10kb_bins.bed"))
  setDT(bin)
  setkey(bin,V1,V2,V3)
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peak <- peak[,c(1:3)]
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  overlaps <- foverlaps(peak, bin, type = "any", nomatch = 0L)  
  
  df <- df[which(df$label %in% overlaps$V4),]
  
  diff <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  diff <- diff[,c("Geneid","Significant")]
  df <- merge(df,diff,by.x="label",by.y="Geneid")
  df <- merge(df,search_table,by.x="sample",by.y="sample_name")
  
  df_avg_percent <- df %>%
    group_by(age, label) %>% 
    summarize(avg_percent = mean(percent, na.rm = TRUE)) 
  df_avg_percent <- as.data.frame(df_avg_percent) 
  young_df <- df_avg_percent[which(df_avg_percent$age=="3M"),]
  old_df <- df_avg_percent[which(df_avg_percent$age=="24M"),]
  colnames(young_df)[3] <- "young"
  colnames(old_df)[3] <- "old"
  t_tissue_summary <- merge(young_df,old_df,by="label")
  t_tissue_summary$delta <- t_tissue_summary$old - t_tissue_summary$young
  t_tissue_summary <- t_tissue_summary[,c("label","young","old","delta")]
  t_tissue_summary <- merge(t_tissue_summary,diff,by.x="label",by.y="Geneid")
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}
to_plot <- tissue_summary
to_plot$delta <- to_plot$delta * 100
to_plot$Significant <- factor(to_plot$Significant,levels=c("Up","Stable","Down"))
p <- ggplot(to_plot, aes(x = Significant, y = delta,fill=Significant)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "DNA methylation", title = "all tissues DNA methylation") +
  theme_bw()+
  ylim(-8,8)
p
ggsave("result/Sup_figures/all_tissues_H3K27me3_bin_young_peak_delta.pdf",p,height = 8,width = 6)

