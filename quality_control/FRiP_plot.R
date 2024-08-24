rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(gg.gap)
# antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me1")
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary","BAT","mammarygland")
FRiP<-data.frame(tissue = character(),  
                  sample_name = character(),
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
      t_FRiP<-data.frame(tissue = tissue,  
                         sample_name = sample_names[k],
                         antibody = antibody,
                         FRiP = df[1,k]/sum(df[,k]),
                         batch = batch[k],  
                         stringsAsFactors = FALSE)  
      FRiP <- rbind(FRiP,t_FRiP)
    }    
  }
}
FRiP <- merge(FRiP,search_table[,c("sample_name","mouse_ID","age")],by="sample_name")
FRiP$FRiP <- FRiP$FRiP*100
FRiP$tissue[which(FRiP$tissue=="brain")] <- "Cortex"
FRiP$tissue[which(FRiP$tissue=="Hip")] <- "Hippocampus"
FRiP$tissue[which(FRiP$tissue=="CB")] <- "Cerebellum"
FRiP <- FRiP %>% mutate(tissue = str_to_title(tissue))  
FRiP$tissue[which(FRiP$tissue=="Bonemarrow")] <- "Bone Marrow"
FRiP$tissue[which(FRiP$tissue=="Bat")] <- "BAT"
FRiP$tissue[which(FRiP$tissue=="Mammarygland")] <- "Mammary Gland"
FRiP$age <- factor(FRiP$age,levels=c("3m","24m"))
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  p<-ggplot(FRiP[which(FRiP$antibody==antibody),],aes(x=tissue,y=FRiP,color = age))+
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    ggtitle(paste0(antibody," FRiP"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("FRiP(%)")
  ggsave(paste0("result/all/QC/FRiP/plot/union/per_antibody/",antibody,"_FRiP_rm_blacklist_age_dotplot.png"),p,width = 15,height = 10,type="cairo")
  # ggsave(paste0("result/all/QC/FRiP/plot/intersect/per_antibody/",antibody,"_FRiP_rm_blacklist_age_dotplot.png"),p,width = 15,height = 10,type="cairo")
  }

for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  p<-ggplot(FRiP[which(FRiP$antibody==antibody),],aes(x=tissue,y=FRiP,color = batch))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    ggtitle(paste0(antibody," FRiP"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("FRiP(%)")
  ggsave(paste0("result/all/QC/FRiP/plot/union/per_antibody/",antibody,"_FRiP_rm_blacklist_batch_dotplot.png"),p,width = 15,height = 10,type="cairo")
  # ggsave(paste0("result/all/QC/FRiP/plot/intersect/per_antibody/",antibody,"_FRiP_rm_blacklist_batch_dotplot.png"),p,width = 15,height = 10,type="cairo")
  }
tissues <- unique(FRiP$tissue)
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  p<-ggplot(FRiP[which(FRiP$tissue==tissue),],aes(x=antibody,y=FRiP,color = age))+
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    ggtitle(paste0(tissue," FRiP"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("FRiP(%)")
  ggsave(paste0("result/all/QC/FRiP/plot/union/per_tissue/",tissue,"_FRiP_rm_blacklist_age_dotplot.png"),p,width = 15,height = 10,type="cairo")
  # ggsave(paste0("result/all/QC/FRiP/plot/intersect/per_tissue/",tissue,"_FRiP_rm_blacklist_age_dotplot.png"),p,width = 15,height = 10,type="cairo")
}

