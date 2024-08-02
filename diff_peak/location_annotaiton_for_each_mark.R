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
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
tissues = c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus")
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
for (i in c(1:length(antibodys))){
  antibody<-antibodys[i]
  up_list <- list()
  down_list <- list()
  for(j in c(1:length(tissues))){
    tissue <- tissues[j]
    out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv"))
    out_up <- out[which(out$Significant_bar=="Up"),]
    out_down <- out[which(out$Significant_bar=="Down"),]
    if(nrow(out_up) > 0){
      up_peak <- GRanges(seqnames = out_up$Chr,   
                           ranges = IRanges(start = out_up$Start, end = out_up$End))
      up_peak_anno <- annotatePeak(up_peak, tssRegion=c(-3000, 3000),
                                     TxDb=txdb, annoDb="org.Mm.eg.db")
      up_list[[j]] <- up_peak_anno
      names(up_list)[j] <- tissue
    }else{
      up_list[[j]] <- NULL
    }
    
    if(nrow(out_down) > 0){
      down_peak <- GRanges(seqnames = out_down$Chr,   
                         ranges = IRanges(start = out_down$Start, end = out_down$End))
      down_peak_anno <- annotatePeak(down_peak, tssRegion=c(-3000, 3000),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")
      down_list[[j]] <- down_peak_anno
      names(down_list)[j] <- tissue
    }else{
      down_list[[j]] <- NULL
    }
  }
  up_list <- Filter(Negate(is.null), up_list)
  plotAnnoBar(up_list)
  ggsave(paste0("result/all/diff/",antibody,"/bins_incease_location_barplot.png"),width = 6,height = 5,type="cairo")
  down_list <- Filter(Negate(is.null), down_list)
  plotAnnoBar(down_list)
  ggsave(paste0("result/all/diff/",antibody,"/bins_decease_location_barplot.png"),width = 6,height = 5,type="cairo")
}
