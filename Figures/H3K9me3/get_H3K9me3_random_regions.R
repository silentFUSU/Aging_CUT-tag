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

## random regions
mouse_genome <- BSgenome.Mmusculus.UCSC.mm10
chrom_lengths <- seqlengths(mouse_genome)
chrom_lengths <- chrom_lengths[which(names(chrom_lengths) %in% paste0("chr",c(c(1:19),"X")))]
min_size <- 200000
regions_size <- 500
selected_regions <- GRanges()
get_random_region <- function(chr_length) {
  start_pos <- sample(1:(chr_length - min_size), 1)
  end_pos <- start_pos + min_size - 1
  return(c(start_pos, end_pos))
}
is_non_overlapping <- function(new_region, selected_regions) {
  overlaps <- findOverlaps(new_region, selected_regions)
  return(length(overlaps) == 0)
}
while (length(selected_regions) < regions_size) {
  for (chr in names(chrom_lengths)) {
    if (length(selected_regions) >= regions_size) break
    
    chr_length <- chrom_lengths[chr]
    
    attempt <- 0
    
    repeat {
      random_region <- get_random_region(chr_length)
      new_range <- GRanges(seqnames = chr,
                           ranges = IRanges(start = random_region[1], end = random_region[2]))
      
      if (is_non_overlapping(new_range, selected_regions)) {
        selected_regions <- c(selected_regions, new_range)
        break
      }
      
      attempt <- attempt + 1
      if (attempt > 100) break  
    }
  }
}
selected_regions <- selected_regions[seq_len(regions_size)]
random_regions <-  data.frame(
  V1 = as.character(seqnames(selected_regions)),
  V2 = start(selected_regions),
  V3 = end(selected_regions)
)
random_regions$V1 <- factor(random_regions$V1,levels=paste0("chr",c(1:19,"X")))
write.table(random_regions[order(random_regions$V1,random_regions$V2),],c("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed"),row.names = F,col.names = F,append = F,sep = "\t",quote = F)

## random regions without overlap with recursion peaks
recursion_peaks <- read.table(c("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed"))
recursion_peaks <- GRanges(seqnames = recursion_peaks$V1,ranges = IRanges(start = recursion_peaks$V2,end = recursion_peaks$V3))
selected_regions <- GRanges()
is_non_overlapping <- function(new_region, selected_regions, recursion_peaks) {
  overlaps1 <- findOverlaps(new_region, selected_regions)
  overlaps2 <- findOverlaps(new_region, recursion_peaks)
  return(length(overlaps1) == 0 && length(overlaps2) == 0)
}
while (length(selected_regions) < regions_size) {
  for (chr in names(chrom_lengths)) {
    if (length(selected_regions) >= regions_size) break
    
    chr_length <- chrom_lengths[chr]
    attempt <- 0
    
    repeat {
      random_region <- get_random_region(chr_length)
      new_range <- GRanges(seqnames = chr,
                           ranges = IRanges(start = random_region[1], end = random_region[2]))
      
      if (is_non_overlapping(new_range, selected_regions, recursion_peaks)) {
        selected_regions <- c(selected_regions, new_range)
        break
      }
      
      attempt <- attempt + 1
      if (attempt > 100) break 
    }
  }
}
selected_regions <- selected_regions[seq_len(regions_size)]
random_regions2 <- data.frame(
  V1 = as.character(seqnames(selected_regions)),
  V2 = start(selected_regions),
  V3 = end(selected_regions)
)
random_regions2$V1 <- factor(random_regions2$V1,levels=paste0("chr",c(1:19,"X")))
write.table(random_regions2[order(random_regions2$V1,random_regions2$V2),],c("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed"),row.names = F,col.names = F,append = F,sep = "\t",quote = F)




