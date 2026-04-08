rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyverse)  
sample <-"WJH-Liver-103"
chr <- "chrY"
tissue <- "liver"
resolution <- "50000"
compartment_recorrect <- function(sample,chr){
  df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
  df$V6[which(df$V2=="chrY")] <- -df$V6[which(df$V2=="chrY")]
  colnames(df) <- c("#peakID","chr","start","end","strand","PC1")
  write.table(df,paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1_recorrect.txt"),row.names = F,quote = F,sep = "\t",append = F)
  
  bedgraph <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.bedGraph"),skip = 1)
  bedgraph$V4[which(bedgraph$V1=="chrY")] <- -bedgraph$V4[which(bedgraph$V1=="chrY")]
  write.table(bedgraph,paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1_recorrect.bedGraph"),row.names = F,quote = F,sep = "\t",append = F,col.names = F)
  
  }