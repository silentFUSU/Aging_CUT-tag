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
bin_size <- "10kb"
tissue_1 <- "spleen"
H3K27me3_1 <- read.csv(paste0("data/samples/",tissue_1,"/H3K27me3/H3K27me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
H3K9me3_1 <- read.csv(paste0("data/samples/",tissue_1,"/H3K9me3/H3K9me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
logFC_1<-merge(H3K27me3_1[which(H3K27me3_1$FDR.old.young < 0.05),c("Geneid","Chr","Start","End","LogFC.old.young")],H3K9me3_1[which(H3K9me3_1$FDR.old.young < 0.05),c("Geneid","LogFC.old.young")],by="Geneid")
colnames(logFC_1)[c(5,6)] <- c(paste0(tissue_1,"_logFC(old-young)_H3K27me3"),paste0(tissue_1,"_logFC(old-young)_H3K9me3"))
logFC_1 <- logFC_1[which(logFC_1[,5]>0 & logFC_1[,6]<0),]

tissue_2 <- "lung"
H3K27me3_2 <- read.csv(paste0("data/samples/",tissue_2,"/H3K27me3/H3K27me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
H3K9me3_2 <- read.csv(paste0("data/samples/",tissue_2,"/H3K9me3/H3K9me3_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
logFC_2<-merge(H3K27me3_2[which(H3K27me3_2$FDR.old.young < 0.05),c("Geneid","Chr","Start","End","LogFC.old.young")],H3K9me3_2[which(H3K9me3_2$FDR.old.young < 0.05),c("Geneid","LogFC.old.young")],by="Geneid")
colnames(logFC_2)[c(5,6)] <- c(paste0(tissue_2,"_logFC(old-young)_H3K27me3"),paste0(tissue_2,"_logFC(old-young)_H3K9me3"))
logFC_2 <- logFC_2[which(logFC_2[,5]>0 & logFC_2[,6]<0),]

logFC <- merge(logFC_1,logFC_2,by="Geneid")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
peak <- GRanges(  
  seqnames = Rle(logFC$Chr.x),  
  ranges = IRanges(start = logFC$Start.x, end = logFC$End.x)
  )  
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno<-unique(as.data.frame(peakAnno))
