rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
library(gridExtra)
library(grid)  
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggsignif)
options(bitmapType="cairo") 
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

tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
tissue_summary <- data.frame()
CpG <- read.table("~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt")
CpG <- CpG[which(CpG$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
CpG <- as.data.table(CpG)
setDT(CpG)
setkey(CpG,V1,V2,V3)
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]

  bin <- as.data.table(read.table("~/ref_data/mm10_10kb_bins.bed"))
  setDT(bin)
  setkey(bin,V1,V2,V3)
  peak <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peak <- peak[,c(1:3)]
  peak <- as.data.table(peak)
  setDT(peak)
  setkey(peak,V1,V2,V3)
  overlaps <- foverlaps(peak, bin, type = "any", nomatch = 0L)  
  
  diff <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  diff <- diff[which(diff$Geneid %in% overlaps$V4),]
  diff_bin <- as.data.table(diff[,c(1:4)])
  setDT(diff_bin)
  setkey(diff_bin,Chr,Start,End)
  
  overlaps <- foverlaps(diff_bin, CpG, type = "any", nomatch = 0L)
  diff$condition <- "out"
  diff$condition[which(diff$Geneid %in% overlaps$Geneid)] <- "CpG island"
  t_tissue_summary <- diff %>%
    group_by(Significant, condition) %>%       
    summarise(count = n(), .groups = "keep") %>%
    mutate(percentage = (count / sum(count)) * 100) 
  t_tissue_summary <- t_tissue_summary %>%
    group_by(Significant) %>%                    
    mutate(percentage = (count / sum(count)) * 100) %>% 
    ungroup() 
  t_tissue_summary <- as.data.frame(t_tissue_summary)
  t_tissue_summary <- t_tissue_summary[which(t_tissue_summary$condition=="CpG island"),]
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
  }
to_plot <- tissue_summary
to_plot$Significant <- factor(to_plot$Significant,levels=c("Up","Stable","Down"))
p <- ggplot(to_plot, aes(x = Significant, y = percentage,fill=Significant)) +
  geom_boxplot(outliers = F) +
  labs(x = NULL, y = "CpG island percentage", title = "CpG island percentage") +
  theme_bw()
ggsave("result/Sup_figures/all_tissues_H3K27me3_bin_CpG_island.pdf",p,width = 6,height = 8)
