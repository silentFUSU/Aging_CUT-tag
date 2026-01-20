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
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$tissue == tissue & search_table$age=="3m" & search_table$antibody=="H3K27me3"),]
  hmr <- read.table(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_young_hmr.counts"),header = T)
  
  blacklist <- read.table("~/ref_data/mm10-blacklist.v2.bed",sep = "\t")
  blacklist <- as.data.table(blacklist)
  setDT(blacklist)
  setkey(blacklist,V1,V2,V3)
  
  hmr_regions <- as.data.table(hmr)
  setDT(hmr_regions)
  setkey(hmr_regions,Chr,Start,End)
  overlaps <- foverlaps(hmr_regions, blacklist, type = "any", nomatch = 0L)  
  hmr <- hmr[which(!hmr$Geneid %in% overlaps$Geneid),]
  
  counts = hmr[,c(7:ncol(hmr))]
  rownames(counts)= hmr$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  counts <- counts[,search_table$sample_name]
  total_reads <- colSums(counts)
  length <- as.numeric(hmr$Length)
  rpkm <- sweep(counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  rpkm$mean <- rowMeans(rpkm)
  rpkm <- rpkm[order(rpkm$mean,decreasing = T),]
  
  
  
  
  top_hmr <- rownames(rpkm)[1:1000]
  regions <- as.data.table(hmr[which(hmr$Geneid %in% top_hmr),c("Geneid","Chr","Start","End")])
  setDT(regions)
  setkey(regions,Chr,Start,End)
  
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  summary <- data.frame()
  for(sample in search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, regions, type = "any", nomatch = 0L)  
    # overlaps <- overlaps[which(overlaps$V5 > 5),]
    
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

write.csv(tissue_summary,"data/samples/WGBS/all_tissues_delta_in_1000_change_in_high_H3K27me3_level_hmr.csv")
