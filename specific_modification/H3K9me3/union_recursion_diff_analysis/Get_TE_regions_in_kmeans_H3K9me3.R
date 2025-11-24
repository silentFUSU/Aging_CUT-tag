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
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggsignif)
options(bitmapType="cairo")  
kmeans <- paste0("kmeans",c(1:4))
TE <- fread("~/ref_data/TE_reference/mm10_TE_metadata_20250116.bed")
setDT(TE)
setkey(TE,V1,V2,V3)
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])

# family_data <-family_data[which(family_data$family_id=="ERV1"),]
# family_data <-family_data[which(family_data$class_id=="LTR"),]
family_data <- as.data.table(family_data)
setDT(family_data)
setkey(family_data,seqnames,start,end)

regions <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/kmeans1_uinon_recursion_peaks.bed"))
regions <- as.data.table(regions)
setDT(regions)
setkey(regions,V1,V2,V3)

overlaps <- foverlaps(family_data, regions, type = "any", nomatch = 0L)  
overlaps <- overlaps[,c("seqnames","start","end")]
write.table(overlaps,"~/ref_data/TE_reference/mm10_TE_kmeans1_regions.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)



write.table(family_data[,c(1:3)],"~/ref_data/TE_reference/mm10_TE_all_regions.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)


