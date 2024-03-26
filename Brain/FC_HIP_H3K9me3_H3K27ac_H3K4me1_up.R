rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
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
library(pheatmap)
library(reshape2)
library(UpSetR)
library(clipr)  
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
tissue <- "brain"
condition1 <- "Down"
condition2 <- "Up"
H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_",bin_size("H3K9me3"),"_bins_diff_after_remove_batch_effect.csv")) 
H3K9me3 <- H3K9me3[which(H3K9me3$Significant_bar==condition1),]
H3K4me1 <- read.csv(paste0("data/samples/",tissue,"/H3K4me1/H3K4me1_",bin_size("H3K4me1"),"_bins_diff_after_remove_batch_effect.csv")) 
H3K4me1 <- H3K4me1[which(H3K4me1$Significant_bar==condition2),]
H3K27ac <- read.csv(paste0("data/samples/",tissue,"/H3K27ac/H3K27ac_",bin_size("H3K27ac"),"_bins_diff_after_remove_batch_effect.csv")) 
H3K27ac <- H3K27ac[which(H3K27ac$Significant_bar==condition2),]
list <- list(H3K9me3_Down=H3K9me3$Geneid,H3K4me1_Up=H3K4me1$Geneid,H3K27ac_Up=H3K27ac$Geneid)
upset(fromList(list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 20, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2 # 文字标签的大小
)
FC <- intersect(H3K9me3$Geneid,H3K4me1$Geneid)
FC <- intersect(FC,H3K27ac$Geneid)
Hip <- intersect(H3K9me3$Geneid,H3K4me1$Geneid)
Hip <- intersect(Hip,H3K27ac$Geneid)
intersect <- intersect(FC,Hip)
FC <- as.data.frame(FC)
Hip <- as.data.frame(Hip)
mm10_10kb <- read.table("/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.bed")
FC <- mm10_10kb[which(mm10_10kb$V4 %in% FC$FC),]
Hip <- mm10_10kb[which(mm10_10kb$V4 %in% Hip$Hip),]

FC_grange_obj <- GRanges(seqnames = FC$V1,   
                         ranges = IRanges(start = FC$V2, end = FC$V3)) 
FC_peakAnno <- annotatePeak(FC_grange_obj, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")


Hip_grange_obj <- GRanges(seqnames = Hip$V1,   
                         ranges = IRanges(start = Hip$V2, end = Hip$V3)) 
Hip_peakAnno <- annotatePeak(Hip_grange_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")

plotAnnoBar(list(FC=FC_peakAnno,Hip=Hip_peakAnno))

FC_peakAnno <- unique(as.data.frame(FC_peakAnno))
Hip_peakAnno <- unique(as.data.frame(Hip_peakAnno))
genelist<-bitr(Hip_peakAnno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <-enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
barplot(genelist_GO,showCategory = 20,font.size=15,label_format = 100)
