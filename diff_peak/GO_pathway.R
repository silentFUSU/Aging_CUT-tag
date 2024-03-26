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
antibody = "H3K4me3"
tissue = "bonemarrow"
# window_size = "1000"
# gap_size = "3000"
# e_value = "100"
# diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
# diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv"))

bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
tissues <- c("spleen","testis","colon","kidney","lung","liver","muscle","Hip","cecum","bonemarrow","brain")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv"))
  diff_up <- diff[which(diff$Significant_bar=="Up"),]
  peak_obj <- GRanges(seqnames = diff_up$Chr,   
                      ranges = IRanges(start = diff_up$Start, end = diff_up$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  # genelist_up <-bitr(diff$symbol[which(diff$Significant_bar=="Up" & str_detect(diff$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  # genelist_up <-bitr(diff$symbol[which(diff$Significant_bar=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  # genelist_up <-bitr(diff$symbol[which(diff$FDR.old.young<0.05 & diff$LogFC.old.young > 0.5)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up <- bitr(peak_anno$SYMBOL[which(str_detect(peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                             OrgDb = GO_database,
                             keyType = "ENTREZID",#设定读取的gene ID类型
                             ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                             pvalueCutoff = 0.05,#设定p值阈值
                             qvalueCutoff = 0.05,#设定q值阈值
                             readable = T)
  if(nrow(genelist_up_GO) > 0){
    barplot(genelist_up_GO,font.size=15,label_format = 100)
    dir.create(paste0("result/",tissue,"/",antibody))
    ggsave(paste0("result/",tissue,"/",antibody,"/promoter_increase_Go.png"),width = 8,height = 7,type="cairo")
  }
  # genelist_down <-bitr(diff$symbol[which(diff$Significant_bar=="Down" & str_detect(diff$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  # genelist_down <-bitr(diff$symbol[which(diff$FDR.old.young<0.05 & diff$LogFC.old.young < -0.5)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  diff_down <- diff[which(diff$Significant_bar=="Down"),]
  peak_obj <- GRanges(seqnames = diff_down$Chr,   
                      ranges = IRanges(start = diff_down$Start, end = diff_down$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  genelist_down <- bitr(peak_anno$SYMBOL[which(str_detect(peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  if(nrow(genelist_down_GO) > 0){
    barplot(genelist_down_GO,font.size=15,label_format = 100)
    dir.create(paste0("result/",tissue,"/",antibody))
    ggsave(paste0("result/",tissue,"/",antibody,"/promoter_decrease_Go.png"),width = 8,height = 7,type="cairo")
  }
}

