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
# antibody <- "ATAC"
antibodys <- c("ATAC","H3K27ac","H3K4me3")
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  tsse <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/tsse/",antibody,"/tsse.txt"), header = FALSE)
  tsse <- tsse[order(tsse$V2),]
  tsse$V1[which(tsse$V1=="brain")] <- "Cortex"
  tsse$V1[which(tsse$V1=="Hip")] <- "Hippocampus"
  tsse$V1[which(tsse$V1=="CB")] <- "Cerebellum"
  tsse <- tsse %>% mutate(V1 = str_to_title(V1))  
  tsse$V1[which(tsse$V1=="Bonemarrow")] <- "Bone Marrow"
  
  tissue <- tissue_antibody(antibody)
  tsse$batch<-rep(c("rep1", "rep1", "rep2", "rep2"), times = length(tissue)) 
  
  ggplot(tsse,aes(x=V1,y=V3,fill = batch))+geom_boxplot()+
    scale_fill_brewer(palette="Set3")+guides(fill = FALSE) +
    ggtitle(paste0(antibody," TSSE"))+ylim(0,15)+
    theme_bw()+theme(text = element_text(size = 18))+xlab("")+labs(fill = "", color = "") +ylab("tsse")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  ggsave(paste0("result/all/QC/tsse/",antibody,"_tsse.png"),width=15,height = 10,type="cairo")
  
}



