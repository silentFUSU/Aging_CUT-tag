rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)

tissues <- c("colon","cecum","ileum","jejunum")
FRiP<-data.frame(sample_name = character(),
                 antibody = character(),
                 FRiP = numeric(),
                 batch = character(),  
                 stringsAsFactors = FALSE)  



batch <- c("batch1","batch1","batch2","batch2")
search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table_ATAC <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(search_table,search_table_ATAC)
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1","ATAC")
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  for (j in c(1:length(tissues))){
    tissue <- tissues[j]
    if(antibody %in% c("H3K27me3","H3K36me3","H3K9me3")){
      df <- read.delim(paste0("result/all/QC/FRiP/union/",tissue,"_",antibody,"_young_old_merge-W1000-G3000-E100_rm_blacklist.counts.summary"),row.names = 1)
      # df <- read.delim(paste0("result/all/QC/FRiP/intersect/",tissue,"_",antibody,"_young_old_intersect-W1000-G3000-E100_rm_blacklist.counts.summary"),row.names = 1)
    }else{
      df <- read.delim(paste0("result/all/QC/FRiP/union/",tissue,"_",antibody,"_macs_young_old_narrowpeak_rm_blacklist.counts.summary"),row.names = 1)
      # df <- read.delim(paste0("result/all/QC/FRiP/intersect/",tissue,"_",antibody,"_macs_young_old_intersect_narrowpeak_rm_blacklist.counts.summary"),row.names = 1)
    }
    if(antibody=="ATAC" & tissue =="testis"){
      df <- df[,-c(1:2)]
    }
    colnames <- colnames(df)
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    sample_names <- gsub(pattern, "\\1", colnames)
    for (k in c(1:ncol(df))){
      t_FRiP<-data.frame(sample_name = sample_names[k],
                         antibody = antibody,
                         FRiP = df[1,k]/sum(df[,k]),
                         batch = batch[k],  
                         stringsAsFactors = FALSE)  
      FRiP <- rbind(FRiP,t_FRiP)
    }    
  }
}

to_plot <- merge(FRiP,search_table[c(1,3,4,5)],by="sample_name")
to_plot <- to_plot[which(to_plot$antibody=="H3K36me3"),]
to_plot <- to_plot[order(to_plot$tissue,to_plot$mouse_ID),]
to_plot$sample_name <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
to_plot$sample_name <- factor(to_plot$sample_name, levels=to_plot$sample_name)
to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
to_plot$FRiP <- to_plot$FRiP*100
ggplot(to_plot,aes(x=sample_name,y=FRiP,color = tissue,shape=age))+
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
  geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
  scale_fill_brewer(palette="Set3")+
  ggtitle(paste0("Intestine H3K36me3 FRiP"))+ylim(0,100)+
  geom_hline(yintercept = mean(to_plot$FRiP), linetype = "dashed", color = "red") +
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("FRiP(%)")
