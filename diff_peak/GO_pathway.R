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
antibody = "H3K36me3"
tissue = "muscle"
window_size = "1000"
gap_size = "3000"
e_value = "100"
diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv"))
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
genelist_up <-bitr(diff$symbol[which(diff$Significant=="Up" & str_detect(diff$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up <-bitr(diff$symbol[which(diff$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up <-bitr(diff$symbol[which(diff$FDR.old.young<0.05 & diff$LogFC.old.young > 0.5)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
barplot(genelist_up_GO,font.size=15,label_format = 100)

genelist_down <-bitr(diff$symbol[which(diff$Significant=="Down" & str_detect(diff$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down <-bitr(diff$symbol[which(diff$FDR.old.young<0.05 & diff$LogFC.old.young < -0.5)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
barplot(genelist_down_GO,font.size=15,label_format = 100)
