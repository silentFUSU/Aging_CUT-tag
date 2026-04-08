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

fai <- read.table("~/ref_data/for_normal_mapping/mm10/mm10.fa.fai", sep="\t", stringsAsFactors=FALSE)
custom.genome <- toGRanges(data.frame(chr=fai$V1, start=1, end=fai$V2))

regions<- read.table("data/samples/brain/H3K27me3/peaks/edd/edd_peaks_fdr05.bed")
granges_object <- GRanges(
  seqnames = regions$V1,
  ranges   = IRanges(start = regions$V2, end = regions$V3)
)
# H3K9me3_peaks <- read.table("data/samples/brain/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100.bed")
# H3K9me3_peaks$length <- H3K9me3_peaks$V3 -H3K9me3_peaks$V2 +1
# H3K9me3_peaks <- H3K9me3_peaks[which(H3K9me3_peaks$length > 100000),]
# peaks_granges_object <- GRanges(
#   seqnames = H3K9me3_peaks$V1,
#   ranges   = IRanges(start = H3K9me3_peaks$V2, end = H3K9me3_peaks$V3)
# )

pdf("result/Sup_figures/H3K27me3_domains_location_mm10.pdf", width = 12, height = 6)  # 可按需要改尺寸
kp <- plotKaryotype(genome = custom.genome, main="H3K27me3 domains", chromosomes = mouse.chromosomes,plot.type=6,cex=1.8)
kpDataBackground(kp, color = "#FFFFFFAA")
kpPlotDensity(kp, granges_object,window.size = 0.5e6, r0=0, r1=0.5,data.panel="ideogram", col="#3c5488", border="#3c5488")
# kpPlotDensity(kp, peaks_granges_object,window.size = 0.5e6, r0=0.51, r1=1,data.panel="ideogram", col="#e64b35", border="#e64b35")

dev.off()
