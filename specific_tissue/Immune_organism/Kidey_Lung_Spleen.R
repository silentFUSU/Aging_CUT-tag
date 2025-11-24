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
tissues <- c("kidney","lung","spleen","bonemarrow","liver")
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
peak_list <- list()
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_",bin_size("H3K9me3"),"_bins_diff_after_remove_batch_effect.csv")) 
  H3K9me3 <- H3K9me3[which(H3K9me3$Significant_bar=="Down"),]
  H3K27me3 <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_",bin_size("H3K27me3"),"_bins_diff_after_remove_batch_effect.csv")) 
  H3K27me3 <- H3K27me3[which(H3K27me3$Significant_bar=="Up"),]
  peak <- intersect(H3K9me3$Geneid, H3K27me3$Geneid)
  peak_list[[i]] <- peak
  names(peak_list)[i] <- tissues[i]
}
upset(fromList(peak_list),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 100,     # 绘制的最大集合个数
      nintersects = 20, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据
      keep.order = T, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE，它根据集合的大小排序。
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 2 # 文字标签的大小
)

mm10_10kb <- read.table("/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.bed")
grange_list <- list()
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- peak_list[[i]]
  df <- mm10_10kb[which(mm10_10kb$V4%in%df),]
  grange_list[[i]] <- GRanges(seqnames = df$V1,   
                              ranges = IRanges(start = df$V2, end = df$V3))
  names(grange_list)[i] <- tissue
}

peakanno_list <- list()
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  peakanno <- annotatePeak(grange_list[[i]], tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakanno_list[[i]] <- peakanno
  names(peakanno_list)[i] <- tissue
}
plotAnnoBar(peakanno_list)

for(i in c(1:length(tissues))){
  peakanno_list[[i]] <- unique(as.data.frame(peakanno_list[[i]]))
}

genelist<-bitr(peakanno_list[[i]]$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <-enrichGO( genelist$ENTREZID,#GO富集分析
                        OrgDb = GO_database,
                        keyType = "ENTREZID",#设定读取的gene ID类型
                        ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                        pvalueCutoff = 0.05,#设定p值阈值
                        qvalueCutoff = 0.05,#设定q值阈值
                        readable = T)
barplot(genelist_GO,showCategory = 10,font.size=15,label_format = 100)
