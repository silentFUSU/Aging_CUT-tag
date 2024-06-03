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
library(UpSetR)
tissues <- c("muscle","thymus","bladder","tongue","ileum")
annotate_list <- list()
mm10_10kb <- read.table("~/ref_data/mm10_10kb_bins.bed")
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
annotate_list_increase <- list()
annotate_list_decrease <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  diff <- read.csv(paste0("data/samples/",tissue,"/H3K36me3/H3K36me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  increase_bar <- unique(diff$Geneid[which(diff$Significant_bar=="Up")])
  decrease_bar <- unique(diff$Geneid[which(diff$Significant_bar=="Down")])
  increase_bar <- mm10_10kb[which(mm10_10kb$V4 %in% increase_bar),]
  decrease_bar <- mm10_10kb[which(mm10_10kb$V4 %in% decrease_bar),]
  increase_bar_grange <- GRanges(seqnames = increase_bar$V1,   
                                          ranges = IRanges(start = increase_bar$V2, end = increase_bar$V3))
  increase_bar_peak_anno <- annotatePeak(increase_bar_grange, tssRegion=c(-3000, 3000),
                                                  TxDb=txdb, annoDb="org.Mm.eg.db")
  
  decrease_bar_grange <- GRanges(seqnames = decrease_bar$V1,   
                                 ranges = IRanges(start = decrease_bar$V2, end = decrease_bar$V3))
  decrease_bar_peak_anno <- annotatePeak(decrease_bar_grange, tssRegion=c(-3000, 3000),
                                         TxDb=txdb, annoDb="org.Mm.eg.db")
  annotate_list_increase[[i]] <- increase_bar_peak_anno 
  annotate_list_decrease[[i]] <- decrease_bar_peak_anno 
  names(annotate_list_increase)[[i]] <- paste0(tissue)
  names(annotate_list_decrease)[[i]] <- paste0(tissue)
  
}
annotate_list_increase[[1]]@annoStat$Feature <- as.vector(annotate_list_increase[[1]]@annoStat$Feature)
annotate_list_increase[[1]]@annoStat$Feature <- factor(annotate_list_increase[[1]]@annoStat$Feature,levels = c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","5' UTR","3' UTR"  ,"1st Exon" ,"Other Exon","1st Intron"  , "Other Intron" ,"Distal Intergenic","Downstream (<=300)" ))
plotAnnoBar(annotate_list_increase,title = "H3K36me3")

annotate_list_decrease[[1]]@annoStat$Feature <- as.vector(annotate_list_decrease[[1]]@annoStat$Feature)
annotate_list_decrease[[1]]@annoStat$Feature <- factor(annotate_list_decrease[[1]]@annoStat$Feature,levels = c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","5' UTR","3' UTR"  ,"1st Exon" ,"Other Exon","1st Intron"  , "Other Intron" ,"Distal Intergenic","Downstream (<=300)" ))
plotAnnoBar(annotate_list_decrease,title = "H3K36me3")
