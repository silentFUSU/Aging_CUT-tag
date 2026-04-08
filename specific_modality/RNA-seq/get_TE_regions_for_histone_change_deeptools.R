rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
library(gridExtra)
library(grid)  
options(bitmapType="cairo")  
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
kmean <- "kmeans1"
df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))

target <- family_data[which(family_data$family_id=="ERV1"),]
# target <- family_data[which(family_data$gene_id=="IAPEy-int"),]
target <- as.data.table(target)
setDT(target)
setkey(target,seqnames,start,end)

head(df)
df <- as.data.table(df)
setDT(df)
setkey(df,V1,V2,V3)
overlaps <- foverlaps(df, target, type = "any", nomatch = 0L)  

write.table(overlaps[,c(1:3)],"~/ref_data/TE_reference/mm10_ERV1_kmeans1_regions.bed",append = F,quote = F,row.names = F,col.names = F,sep = "\t")
