rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(karyoploteR)
library(GenomicRanges)  
mouse.chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                       "chr6", "chr7", "chr8", "chr9", "chr10",
                       "chr11", "chr12", "chr13", "chr14", "chr15",
                       "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
kmean <- "kmeans4"
df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
df$label <- kmean
granges_object <- GRanges(
  seqnames = df$V1,
  ranges   = IRanges(start = df$V2, end = df$V3),
  label    = df$label
)
kp <- plotKaryotype(genome = "mm10", main=paste0(kmean," H3K9me3 peaks"), chromosomes = mouse.chromosomes,plot.type=6,cex=1.8)
kpDataBackground(kp, color = "#FFFFFFAA")
kpPlotDensity(kp, granges_object,window.size = 0.5e6, data.panel="ideogram", col="#3388FF", border="#3388FF")

