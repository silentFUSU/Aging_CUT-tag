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
library(scales)
library(colorspace)   
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
cecum <- read.delim("data/samples/cecum/H3K27me3_H3K9me3_intersect/10kb_all_significant_fourth_quadrant.bed",header = F)
colon <- read.delim("data/samples/colon/H3K27me3_H3K9me3_intersect/10kb_all_significant_fourth_quadrant.bed",header = F)
intersect <- intersect(cecum$V4,colon$V4)
intersect <- cecum[which(cecum$V4 %in% intersect),]
grange_obj <- GRanges(seqnames = intersect$V1,   
              ranges = IRanges(start = intersect$V2, end = intersect$V3)) 
peakAnno <- annotatePeak(grange_obj, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoBar(peakAnno)
peakAnno<-unique(as.data.frame(peakAnno))
genelist<-bitr(peakAnno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <-enrichGO( genelist$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
barplot(genelist_GO,showCategory = 20,font.size=15,label_format = 100)

grange_obj <- GRanges(seqnames = cecum$V1,   
                      ranges = IRanges(start = cecum$V2, end = cecum$V3)) 
peakAnno <- annotatePeak(grange_obj, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoBar(peakAnno)


grange_obj <- GRanges(seqnames = colon$V1,   
                      ranges = IRanges(start = colon$V2, end = colon$V3)) 
peakAnno <- annotatePeak(grange_obj, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoBar(peakAnno)

