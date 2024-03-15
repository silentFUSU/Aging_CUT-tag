rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(DiffBind)
library(tidyverse)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
sample <- read.table("data/samples/brain/H3K27me3/SampleSheet.csv",header = TRUE, sep = ",", quote = "\"")
dbObj <- dba(sampleSheet=sample)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
dba.plotPCA(dbObj, contrast=1, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, method=DBA_EDGER)
res_deseq <- as.data.frame(dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1))


res_deseq$Significant <- "Stable"
res_deseq$Significant[which(res_deseq$Fold < 0 & res_deseq$FDR < 0.05,)] <- "Down"
res_deseq$Significant[which(res_deseq$Fold > 0 & res_deseq$FDR < 0.05,)] <- "Up"

colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
res_deseq$log2FC <- log2(abs(res_deseq$Fold))

ggplot(
  # 数据、映射、颜色
  res_deseq, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour[[nrow(as.data.frame(table(res_deseq$Significant)))]]) +
  # scale_color_manual(values = c("blue","grey")) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)") +
  # xlim(-3,3)+
  # 图例
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()+theme(text = element_text(size = 18))
