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
search_table_cut <- read.csv("data/samples/all/CUTTag_search_table.csv")
search_table_atac <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(search_table_cut,search_table_atac)
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  tsse <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/tsse/",antibody,"/tsse.txt"), header = FALSE)
  pattern <- ".*(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  tsse$V2 <- gsub(pattern, "\\1", tsse$V2)
  colnames(tsse)[2] <- "sample_name"
  tsse <- merge(search_table,tsse[,c(2:3)],by="sample_name")
 
  tsse$age <- factor(tsse$age, levels = c("3m","24m"))
  ggplot(tsse,aes(x=tissue,y=V3,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = mouse_ID), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_fill_brewer(palette="Set3")+
    geom_hline(yintercept = 4, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 5, color = "blue", linetype = "dashed") + 
    annotate("text", x = Inf, y = 3.5, label = "tsse=4", color = "red", vjust = -0.5, hjust = 1.1, size = 5) +   
    annotate("text", x = Inf, y = 4.5, label = "tsse=5", color = "blue", vjust = -0.5, hjust = 1.1, size = 5) + 
    ggtitle(paste0(antibody," TSSe"))+ylim(0,20)+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("TSSe")
  
  ggsave(paste0("result/all/QC/tsse/",antibody,"_tsse.png"),width=15,height = 10,type="cairo")
  
}



