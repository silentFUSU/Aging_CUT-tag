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
  df <- df[which(df$total_V5 > 15),]
  t_df <- df[which(df$sample %in% search_table$sample_name[which(search_table$age=="3M")]),]
  average_percent_by_label <- t_df %>%
    group_by(label) %>%
    summarise(mean_percent = mean(percent))
  average_percent_by_label <- average_percent_by_label[order(average_percent_by_label$mean_percent),]
  bin <- as.data.table(read.table("~/ref_data/mm10_10kb_bins.bed"))
  setDT(bin)
  setkey(bin,V1,V2,V3)
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peak <- peak[,c(1:3)]
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  overlaps <- foverlaps(peak, bin, type = "any", nomatch = 0L)  
  average_percent_by_label <- as.data.frame(average_percent_by_label[which(average_percent_by_label$label %in% overlaps$V4),])
  average_percent_by_label <- average_percent_by_label[order(average_percent_by_label$mean_percent),]
  LMR <- average_percent_by_label$label[1:1000]
  # LMR <- average_percent_by_label$label[which(average_percent_by_label$mean_percent < 0.2)]
  
  LMR_meth <- df[which(df$label %in% LMR),]
  total_sums_by_tissue_sample <- LMR_meth %>%
    group_by(tissue, sample) %>%
    summarise(
      total_V4_sum = sum(total_V4),
      total_V5_sum = sum(total_V5)
    )
  total_sums_by_tissue_sample$methylation <- total_sums_by_tissue_sample$total_V4_sum/total_sums_by_tissue_sample$total_V5_sum *100
  t_tissue_summary <- data.frame(tissue=tissue_label_change(tissue),
             young_methylation=mean(total_sums_by_tissue_sample$methylation[which(total_sums_by_tissue_sample$sample %in% search_table$sample_name[which(search_table$age=="3M")])]),
             old_methylation=mean(total_sums_by_tissue_sample$methylation[which(total_sums_by_tissue_sample$sample %in% search_table$sample_name[which(search_table$age=="24M")])]),
             delta=mean(total_sums_by_tissue_sample$methylation[which(total_sums_by_tissue_sample$sample %in% search_table$sample_name[which(search_table$age=="24M")])])-mean(total_sums_by_tissue_sample$methylation[which(total_sums_by_tissue_sample$sample %in% search_table$sample_name[which(search_table$age=="3M")])]))
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
  }

tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
to_plot <- tissue_summary
to_plot$tissue <- factor(to_plot$tissue,levels = tissue_order)
to_plot$condition <- "Up"
to_plot$condition[which(to_plot$delta < 0 )] <- "Down"
color <- setNames(c("#f39b7f","#4dbbd5"),c("Up","Down"))
ggplot(to_plot, aes(x = tissue, y = delta, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Delta")+
  ylim(-18,18) +
  ggtitle(NULL) +
  guides(fill = FALSE) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")

write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_1000_LMR_in_H3K27me3_young_peaks.csv")
