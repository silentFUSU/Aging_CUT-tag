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
library(tidyr)
library(MASS) 
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  dir.create(paste0("result/all/diff/",antibody))
  diff <- read.csv(paste0("data/samples/all/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_use_tissue_as_batch_effect.csv"))
  percent <- as.data.frame(table(diff$Significant))
  pie(percent$Freq,labels = percent$Var1)
  ggsave(paste0("result/all/diff/",antibody,"/",antibody,"_pie_percent.png"),width = 10,height = 10,type="cairo")
  diff_down <- diff[which(diff$Significant=="Down"),]
  diff_up <- diff[which(diff$Significant=="Up"),]
  down_peak <- GRanges(seqnames = diff_down$Chr,   
                        ranges = IRanges(start = diff_down$Start, end = diff_down$End))
  down_peak_anno <- annotatePeak(down_peak, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  
  up_peak <- GRanges(seqnames = diff_up$Chr,   
                       ranges = IRanges(start = diff_up$Start, end = diff_up$End))
  up_peak_anno <- annotatePeak(up_peak, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
  plotAnnoBar(list(down=down_peak_anno,up=up_peak_anno))
  
  down_peak_anno <- unique(as.data.frame(down_peak_anno))
  genelist_down <- bitr(down_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  if(nrow(genelist_down_GO) > 0){
    barplot(genelist_down_GO,font.size=15,label_format = 100)
    ggsave(paste0("result/all/diff/",antibody,"/all_tissue_batch_effect",bin_size(antibody),"_decrease_Go.png"),width = 8,height = 7,type="cairo")
  }
  up_peak_anno <- unique(as.data.frame(up_peak_anno))
  genelist_up <- bitr(up_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  if(nrow(genelist_up_GO) > 0){
    barplot(genelist_up_GO,font.size=15,label_format = 100)
    ggsave(paste0("result/all/diff/",antibody,"/all_tissue_batch_effect",bin_size(antibody),"_increase_Go.png"),width = 8,height = 7,type="cairo")
  }
  
}
