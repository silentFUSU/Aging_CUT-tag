rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(ChIPseeker)
fourth <- read.csv("data/samples/MEF_OE/H3K27me3/correlation_with_senescence_fourth_quadrant_common_genes.csv")
genes <- fourth$gene[which(fourth$n==3)]

conditions <- c("Bmi1","Cbx2","Cbx7","mCbx8")
summary <- data.frame()
for(condition in conditions){
  df <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_MEF_",condition,"_diff_expression_gene_strict_filter_bar.csv"))
  df <- df[which(df$Geneid %in% genes),c("Geneid","LogFC.oe.vec","FDR.oe.vec","Significant")]
  colnames(df) <- c("Geneid",paste0(condition," logFC"),paste0(condition," FDR"),paste0(condition," Significant"))
  
  if(nrow(summary)==0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="Geneid",all=T)
  }
}
df <- read.csv(paste0("data/samples/RNA/MEF_mid_age/diff_expression_gene_strict_filter_bar.csv"))
df <- df[which(df$X %in% genes),c("X","logFC","fdr","Significant")]
condition <- "senescence"
colnames(df) <- c("Geneid",paste0(condition," logFC"),paste0(condition," FDR"),paste0(condition," Significant"))
summary <- merge(summary,df,by="Geneid",all=T)
write.csv(summary,"data/samples/RNA/MEF_OE_RNA/four_condition_gene_change_in_fourth_quadrant_genes.csv")
