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
library(patchwork)
bin_size <-"1kb"
antibody <- "ATAC"
peak_preprocess_bin_level <- function(tissue,antibody,bin_size){
  tab = read.delim(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames <- colnames(tab)[7:length(tab)]
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[7:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames] 
  if(tissue=="ovary"){
    age <- c("old","young","old","young")
    colnames(counts) <- c("old_1","young_1","old_2","young_2")
  }else if(tissue=="brain"){
    age <- c("young","young","old","old")
    colnames(counts) <- c("young_1","young_2","old_1","old_2")
  }else{
    age <- c("young","old","young","old")
    colnames(counts) <- c("young_1","old_1","young_2","old_2")
  }
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
  y$samples$group <- factor(y$samples$group,c("young","old"))
  y <- calcNormFactors(y)
  batch <- factor(y$samples$batch)
  if(tissue=="brain"){
    design <- model.matrix(~group, y$samples)
  }else{
    design <- model.matrix(~batch+group, y$samples)
  }


  
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  if(tissue=="brain"){
    lrt = glmLRT(fit_tag, coef = 2)
  }else{
    lrt = glmLRT(fit_tag, coef = 3)
  }

  tab<-tab[keep,]
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  out$Significant_bar <- "Stable"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 > 1.2) & (out$old_2/out$young_2 > 1.2))] <- "Up"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 < 0.8) & (out$old_2/out$young_2 < 0.8))] <- "Down"
  
  write.csv(out,paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"),row.names = F)
  outup <- out[which(out$Significant_bar=="Up"),]
  outdown <- out[which(out$Significant_bar=="Down"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

antibodys <- c("ATAC")
# tissues <- c("stomach","skin")
# tissues <- c("aorta","tongue")
tissues <- c("bladder")
bin_size <-"10kb"
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    peak_preprocess_bin_level(tissue,antibody,bin_size)
  }
}
bin_size <-"1kb"
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    peak_preprocess_bin_level(tissue,antibody,bin_size)
  }
}
