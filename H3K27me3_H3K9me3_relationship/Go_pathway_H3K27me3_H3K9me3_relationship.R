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
antibody = "H3K27me3"
tissue = "brain"
window_size = "1000"
gap_size = "3000"
e_value = "100"
if(antibody=="H3K27me3"){
  another_antibody <- "H3K9me3"
}else{
  another_antibody <- "H3K27me3"
}
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_diff_after_remove_batch_effect_relationship_with_",another_antibody,".csv"))
genelist <- bitr(diff$symbol[which(diff$Significant=="Up" & diff$relationship_between_H3K27me3_H3K9me3==paste0("decrease in ",another_antibody))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <-enrichGO(genelist$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
barplot(genelist_GO,font.size=15,label_format = 100)


