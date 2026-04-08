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
  print(tissue)
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
  
  EZH2 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/EZH2_peaks.narrowPeak")
  EZH2 <- EZH2[,c(1:3,9)]
  colnames(EZH2)[4] <- "EZH2"
  
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
  
  regions <- as.data.table(hmr[which(hmr$V4 %in% score$V4[1:min(1000,nrow(score))]),c(1:4)])
  
 
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody=="H3K27me3"),]
  
  tab <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_all_samples_hmr.counts"),header = T)
  tab_summary <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_all_samples_hmr.counts.summary"),header = T)
  
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  tab_summary <- tab_summary[-2,search_table$sample_name]
  total_reads <- colSums(tab_summary)
  length <- as.numeric(tab$Length)
  rpkm <- sweep(counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  rpkm_young <- rpkm[,search_table$sample_name[which(search_table$age=="3m")]]
  rpkm_old <- rpkm[,search_table$sample_name[which(search_table$age=="24m")]]
  rpkm_young$mean_young <- rowMeans(rpkm_young)
  rpkm_old$mean_old <- rowMeans(rpkm_old)
  rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
  rpkm_mean_summary <- rpkm_mean_summary[which(rpkm_mean_summary$Row.names %in% regions$V4),]
  rpkm_mean_summary$tissue <- tissue_label_change(tissue)
  
  tissue_summary <- rbind(tissue_summary,rpkm_mean_summary)
}
tissue_order <- c("Thymus","Muscle","Cerebellum","Testis","Cortex","Bladder","Hippocampus","Heart","Liver",
                  "Bone Marrow","Lung","Kidney","Aorta","iWAT","BAT","Spleen","Mammary Gland","Stomach","Ovary",
                  "Skin","Uterus","Tongue","Pancreas","Ileum","Jejunum","Cecum","Colon")
to_plot <- reshape2::melt(tissue_summary)
to_plot$tissue <- factor(to_plot$tissue,levels=tissue_order)
p <- ggplot(to_plot, aes(x = tissue, y = value, fill = variable)) +  
  geom_boxplot(outliers = F) +  
  theme_bw() +   
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("H3K27me3 RPKM")+
  ggtitle(NULL) 
p

to_plot$tissue_condition <- "medium"
to_plot$tissue_condition[which(to_plot$tissue %in% c("Thymus","Muscle","Cerebellum","Testis"))] <- "bottom"
to_plot$tissue_condition[which(to_plot$tissue %in% c("Pancreas","Ileum","Jejunum","Cecum","Colon"))] <- "top"
p <- ggplot(to_plot, aes(x = tissue_condition, y = value, fill = variable)) +  
  geom_boxplot(outliers = F) +  
  theme_bw() +   
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("H3K27me3 RPKM")+
  ggtitle(NULL) 
p
