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
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow")
FRiP<-data.frame(tissue = character(),  
                  antibody = character(),
                  FRiP = numeric(),
                  age = character(),
                  batch = character(),  
                  stringsAsFactors = FALSE)  
age<-c("young","old","young","old")
batch <- c("batch1","batch1","batch2","batch2")
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  for (j in c(1:length(tissues))){
    tissue <- tissues[j]
    df <- read.delim(paste0("result/all/QC/FRiP/",tissue,"_",antibody,"_rm_blacklist.counts.summary"),row.names = 1)
    for (k in c(1:ncol(df))){
      t_FRiP<-data.frame(tissue = tissue,  
                         antibody = antibody,
                         FRiP = df[1,k]/sum(df[,k]),
                         age = age[k],
                         batch = batch[k],  
                         stringsAsFactors = FALSE)  
      FRiP <- rbind(FRiP,t_FRiP)
    }    
    }
}

FRiP$FRiP <- FRiP$FRiP*100
FRiP$tissue[which(FRiP$tissue=="brain")] <- "brain_FC"
FRiP$tissue[which(FRiP$tissue=="Hip")] <- "brain_HIP"
FRiP$tissue <- factor(FRiP$tissue,levels=c("brain_FC","brain_HIP","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow"))
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  p<-ggplot(FRiP[which(FRiP$antibody==antibody),],aes(x=tissue,y=FRiP,fill = batch))+geom_boxplot()+
    scale_fill_brewer(palette="Set3")+
    ggtitle(paste0(antibody," FRiP"))+ylim(0,100)+
    theme_bw()+theme(text = element_text(size = 18))+xlab("")+labs(fill = "", color = "") +ylab("FRiP(%)")
  ggsave(paste0("result/all/QC/FRiP/",antibody,"_FRiP_rm_blacklist.png"),p,width = 15,height = 10)
  }
