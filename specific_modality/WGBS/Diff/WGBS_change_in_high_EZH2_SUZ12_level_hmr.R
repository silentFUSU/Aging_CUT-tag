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
tissue_overlap_regions <- data.frame()
for(tissue in tissues){
  print(tissue)
  # hmr <- read.table(paste0("data/samples/WGBS/",tissue,"/hmr/young_hmr.bed"))
  hmr <- read.table(paste0("data/samples/WGBS/",tissue,"/hmr/all_samples_hmr.bed"))
  if(tissue %in% c("mammarygland","ovary","uterus")){
    hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X"))),]
  }else{
    hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  }

  blacklist <- read.table("~/ref_data/mm10-blacklist.v2.bed",sep = "\t")
  blacklist <- as.data.table(blacklist)
  setDT(blacklist)
  setkey(blacklist,V1,V2,V3)
  
  hmr_regions <- as.data.table(hmr)
  setDT(hmr_regions)
  setkey(hmr_regions,V1,V2,V3)
  overlaps <- foverlaps(hmr_regions, blacklist, type = "any", nomatch = 0L)  
  
  hmr <- hmr[which(!hmr$V4 %in% overlaps$i.V4),]
  
  # EZH2 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak/EZH2_peaks.narrowPeak")
  # EZH2 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_p_05/EZH2_peaks.narrowPeak")
  EZH2 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/EZH2_peaks.narrowPeak")
  EZH2 <- EZH2[,c(1:3,9)]
  colnames(EZH2)[4] <- "EZH2"
  
  # SUZ12 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak/SUZ12_peaks.narrowPeak")
  # SUZ12 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_p_05/SUZ12_peaks.narrowPeak")
  SUZ12 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/SUZ12_peaks.narrowPeak")
  SUZ12 <- SUZ12[,c(1:3,9)]
  colnames(SUZ12)[4] <- "SUZ12"
  
  hmr_regions <- as.data.table(hmr)
  setDT(hmr_regions)
  setkey(hmr_regions,V1,V2,V3)
  
  EZH2 <- as.data.table(EZH2)
  setDT(EZH2)
  setkey(EZH2,V1,V2,V3)
  
  SUZ12 <- as.data.table(SUZ12)
  setDT(SUZ12)
  setkey(SUZ12,V1,V2,V3)
  
  overlaps <- foverlaps(EZH2,hmr_regions, type = "any", nomatch = 0L)  
  overlaps <- as.data.table(overlaps[,c("V1","V2","V3","V4","EZH2")])
  setDT(overlaps)
  setkey(overlaps,V1,V2,V3)
  
  overlaps_SUZ12 <- foverlaps(overlaps, SUZ12, type = "any", nomatch = 0L)  
  overlaps_SUZ12$SUZ12 <- as.numeric(overlaps_SUZ12$SUZ12)
  overlaps_SUZ12$EZH2 <- as.numeric(overlaps_SUZ12$EZH2)
  overlaps_SUZ12$score <- (overlaps_SUZ12$SUZ12 + overlaps_SUZ12$EZH2)/2
  score <- overlaps_SUZ12 %>%
    group_by(V4) %>%
    summarize(
      score = max(score, na.rm = TRUE)
    )
  score <- as.data.frame(score)
  score <- score[order(score$score,decreasing = T),]
  t_tissue_overlap_regions <- data.frame(tissue=tissue_label_change(tissue),counts=nrow(score))
  tissue_overlap_regions <- rbind(tissue_overlap_regions,t_tissue_overlap_regions)
  
  score_top <- score[1:min(1000,nrow(score)),]
  score_top <- merge(score_top,hmr[,c(1:4)],by="V4")
  score_top <- score_top[,c("V1","V2","V3","V4","score")]
  score_top <- score_top[order(score_top$score,decreasing = T),]
  write.table(score_top,paste0("data/samples/WGBS/",tissue,"/hmr/PRC2_LMR_top1000.bed"),append = F,quote = F,row.names = F,col.names = F)
  regions <- as.data.table(hmr[which(hmr$V4 %in% score$V4[1:min(1000,nrow(score))]),c(1:3)])
  # regions <- as.data.table(overlaps_SUZ12[,c(1:3)])
  setDT(regions)
  setkey(regions,V1,V2,V3)
  
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  summary <- data.frame()
  for(sample in search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, regions, type = "any", nomatch = 0L)  
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5))]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum*100
    result <- result[,c("V4_sum","V5_sum","methylation")]
    result$sample <- sample 
    summary <- rbind(summary,result)
  }
  t_tissue_summary <- data.frame(tissue=tissue_label_change(tissue),
                                 young_methylation=mean(summary$methylation[which(summary$sample %in% search_table$sample_name[which(search_table$age=="3M")])]),
                                 old_methylation=mean(summary$methylation[which(summary$sample %in% search_table$sample_name[which(search_table$age=="24M")])]),
                                 delta=mean(summary$methylation[which(summary$sample %in% search_table$sample_name[which(search_table$age=="24M")])])-mean(summary$methylation[which(summary$sample %in% search_table$sample_name[which(search_table$age=="3M")])]))
  
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}
tissue_summary <- read.csv("data/samples/WGBS/all_tissues_all_samples_delta_in_top1000_hmr_change_in_high_EZH2_SUZ12_peak_q01_input_level_hmr.csv")

tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
to_plot <- tissue_summary
# to_plot$tissue <- factor(to_plot$tissue,levels = tissue_order)
to_plot <- to_plot[order(to_plot$delta,decreasing = T),]
to_plot$tissue <- factor(to_plot$tissue,levels = rev(to_plot$tissue))
to_plot$condition <- "Up"
to_plot$condition[which(to_plot$delta < 0 )] <- "Down"
color <- setNames(c("#f39b7f","#4dbbd5"),c("Up","Down"))
p <- ggplot(to_plot, aes(x = tissue, y = delta, fill = condition)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Delta")+
  ylim(-1,16) +
  ggtitle(NULL) +
  guides(fill = FALSE) +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")
ggsave("result/figures/WGBS_delta_barplot_in_PRC2.pdf",p,width = 6,height = 4)
# write.csv(tissue_summary,"data/samples/WGBS/all_tissues_all_samples_delta_in_top1000_hmr_change_in_high_EZH2_SUZ12_peak_q01_input_level_hmr.csv")
# to_plot <- tissue_overlap_regions
# to_plot$tissue <- factor(to_plot$tissue,levels = tissue_order)
# 
# ggplot(to_plot, aes(x = tissue, y = counts)) +  
#   geom_bar(stat = 'identity') +   
#   theme_bw() +   
#   scale_fill_manual(values = color) +
#   theme(axis.title.x = element_blank(), 
#         axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
#         text = element_text(size = 20),legend.title = element_blank()) +
#   ylab("counts")+
#   ggtitle(NULL) +
#   guides(fill = FALSE) +
#   geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dashed")

