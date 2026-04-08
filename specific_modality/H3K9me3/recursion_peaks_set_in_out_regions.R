rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
options(scipen = 999)
peaks <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.bed")
chrom_size <- read.table("~/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes")
chrom_size <- chrom_size[which(chrom_size$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
chrom_size$V3 <- chrom_size$V2
chrom_size$V2 <- 1
get_uncovered_regions <- function(chrom, start, end, peaks) {
  chrom_peaks <- peaks %>% filter(V1 == chrom) %>% arrange(V2)
  uncovered <- list()
  if (start < chrom_peaks$V2[1]) {
    uncovered <- append(uncovered, list(c(start, chrom_peaks$V2[1] - 1)))
  }
  for (i in 1:(nrow(chrom_peaks) - 1)) {
    if (chrom_peaks$V3[i] < chrom_peaks$V2[i + 1]) {
      uncovered <- append(uncovered, list(c(chrom_peaks$V3[i] + 1, chrom_peaks$V2[i + 1] - 1)))
    }
  }
  if (chrom_peaks$V3[nrow(chrom_peaks)] < end) {
    uncovered <- append(uncovered, list(c(chrom_peaks$V3[nrow(chrom_peaks)] + 1, end)))
  }
  return(do.call(rbind, uncovered))
}
uncovered_regions <- data.frame(V1 = character(), V2 = numeric(), V3 = numeric())
for (i in 1:nrow(chrom_size)) {
  chrom <- chrom_size$V1[i]
  start <- chrom_size$V2[i]
  end <- chrom_size$V3[i]
  
  uncovered <- get_uncovered_regions(chrom, start, end, peaks)
  if (!is.null(uncovered) && nrow(uncovered) > 0) {
    uncovered_regions <- rbind(uncovered_regions, data.frame(V1 = chrom, V2 = uncovered[, 1], V3 = uncovered[, 2]))
  }
}
peaks$V4 <- paste0("peak",c(1:nrow(peaks)))
uncovered_regions$V4 <- paste0("outpeak",c(1:nrow(uncovered_regions)))
bed <- rbind(peaks,uncovered_regions)
bed$V1 <- factor(bed$V1,levels = paste0("chr",c(1:19,"X","Y")))
bed <- bed %>%
  arrange(V1, V2, V3)
write.table(bed,"data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_all_regions.bed",sep = "\t",col.names = F,row.names = F,append = F,quote = F)
