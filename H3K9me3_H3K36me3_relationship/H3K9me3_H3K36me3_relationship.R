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
tissues <- c("muscle","liver","spleen","colon","cecum")
mm10_10kb <- read.table("~/ref_data/mm10_10kb_bins.bed")
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
annotate_list <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  H3K9me3 <-  read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K36me3 <- read.csv(paste0("data/samples/",tissue,"/H3K36me3/H3K36me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_increase_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Up")])
  H3K9me3_decrease_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Down")])
  H3K36me3_increase_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Up")])
  H3K36me3_decrease_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Down")])
  
  H3K9me3_H3K36me3_down <- intersect(H3K9me3_decrease_bar,H3K36me3_decrease_bar)
  H3K9me3_H3K36me3_up <- intersect(H3K9me3_increase_bar,H3K36me3_increase_bar)
  
  H3K9me3_H3K36me3_down <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_down),]
  H3K9me3_H3K36me3_up <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_up),]
  
  if(nrow(H3K9me3_H3K36me3_down)>0){
    H3K9me3_H3K36me3_down_grange <- GRanges(seqnames = H3K9me3_H3K36me3_down$V1,   
                                            ranges = IRanges(start = H3K9me3_H3K36me3_down$V2, end = H3K9me3_H3K36me3_down$V3))
    H3K9me3_H3K36me3_down_peak_anno <- annotatePeak(H3K9me3_H3K36me3_down_grange, tssRegion=c(-3000, 3000),
                                                    TxDb=txdb, annoDb="org.Mm.eg.db")
  }else{
    H3K9me3_H3K36me3_down_peak_anno <- NULL
  }

  if(nrow(H3K9me3_H3K36me3_up)>0){
  H3K9me3_H3K36me3_up_grange <- GRanges(seqnames = H3K9me3_H3K36me3_up$V1,   
                                         ranges = IRanges(start = H3K9me3_H3K36me3_up$V2, end = H3K9me3_H3K36me3_up$V3))
  H3K9me3_H3K36me3_up_peak_anno <- annotatePeak(H3K9me3_H3K36me3_up_grange, tssRegion=c(-3000, 3000),
                                                 TxDb=txdb, annoDb="org.Mm.eg.db")
  }else{
    H3K9me3_H3K36me3_up_peak_anno <- NULL
  }
  annotate_list[[i*2-1]] <- H3K9me3_H3K36me3_up_peak_anno
  annotate_list[[i*2]] <- H3K9me3_H3K36me3_down_peak_anno 
  names(annotate_list)[[i*2-1]] <- paste0(tissue,"_up")
  names(annotate_list)[[i*2]] <- paste0(tissue,"_down")

}
annotate_list <- Filter(Negate(is.null), annotate_list)
plotAnnoBar(annotate_list,title = "H3K9me3 and H3K36me3")
number_annotate_list <- vector()
for(i in c(1:length(annotate_list))){
  number_annotate_list[i] <- annotate_list[[i]]@peakNum
}

tissues <- c("muscle","liver","spleen")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  H3K9me3 <-  read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K36me3 <- read.csv(paste0("data/samples/",tissue,"/H3K36me3/H3K36me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_increase_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Up")])
  H3K9me3_decrease_bar <- unique(H3K9me3$Geneid[which(H3K9me3$Significant_bar=="Down")])
  H3K36me3_increase_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Up")])
  H3K36me3_decrease_bar <- unique(H3K36me3$Geneid[which(H3K36me3$Significant_bar=="Down")])
  
  H3K9me3_H3K36me3_down <- intersect(H3K9me3_decrease_bar,H3K36me3_decrease_bar)
  H3K9me3_H3K36me3_up <- intersect(H3K9me3_increase_bar,H3K36me3_increase_bar)
  
  H3K9me3_H3K36me3_down <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_down),]
  H3K9me3_H3K36me3_up <- mm10_10kb[which(mm10_10kb$V4 %in% H3K9me3_H3K36me3_up),]
  H3K9me3_H3K36me3_down_grange <- GRanges(seqnames = H3K9me3_H3K36me3_down$V1,   
                                          ranges = IRanges(start = H3K9me3_H3K36me3_down$V2, end = H3K9me3_H3K36me3_down$V3))
  H3K9me3_H3K36me3_down_peak_anno <- annotatePeak(H3K9me3_H3K36me3_down_grange, tssRegion=c(-3000, 3000),
                                                  TxDb=txdb, annoDb="org.Mm.eg.db")
  H3K9me3_H3K36me3_up_grange <- GRanges(seqnames = H3K9me3_H3K36me3_up$V1,   
                                        ranges = IRanges(start = H3K9me3_H3K36me3_up$V2, end = H3K9me3_H3K36me3_up$V3))
  H3K9me3_H3K36me3_up_peak_anno <- annotatePeak(H3K9me3_H3K36me3_up_grange, tssRegion=c(-3000, 3000),
                                                TxDb=txdb, annoDb="org.Mm.eg.db")
  
  H3K9me3_H3K36me3_up_peak_anno <- unique(as.data.frame(H3K9me3_H3K36me3_up_peak_anno))

  genelist_up <- bitr(H3K9me3_H3K36me3_up_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  barplot(genelist_up_GO,font.size=15,label_format = 100)
  H3K9me3_H3K36me3_down_peak_anno <- unique(as.data.frame(H3K9me3_H3K36me3_down_peak_anno))
  
  genelist_down <- bitr(H3K9me3_H3K36me3_down_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  barplot(genelist_down_GO,font.size=15,label_format = 100)
}