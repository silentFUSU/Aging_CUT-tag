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
regions<- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_chr <- strsplit(as.character(regions$label), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  

regions <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=regions$cluster)
regions$start <- as.numeric(regions$start)
regions$end <- as.numeric(regions$end)

kmean <- "Stable"
df <- regions[which(regions$cluster==kmean),]
granges_object <- GRanges(
  seqnames = df$chr,
  ranges   = IRanges(start = df$start, end = df$end),
  label    = df$cluster
)

kp <- plotKaryotype(genome = "mm10", main=paste0(kmean," H3K9me3 peaks"), chromosomes = mouse.chromosomes,plot.type=6,cex=1.8)
kpDataBackground(kp, color = "#FFFFFFAA")
kpPlotDensity(kp, granges_object,window.size = 0.5e6, data.panel="ideogram", col="#3388FF", border="#3388FF")

# regions <- createRandomRegions(nregions=400, length.mean = 3e6, mask=NA, non.overlapping = FALSE)
kp <- plotKaryotype(genome = "mm10", main=paste0(kmean," H3K9me3 peaks"), chromosomes = mouse.chromosomes,cex=1.8)
kpPlotRegions(kp, data=granges_object,avoid.overlapping = T)
