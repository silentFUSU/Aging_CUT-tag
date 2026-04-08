rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(GenomicRanges)

tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}

read_bed_to_granges <- function(file_path) {
  bed_data <- read.table(file_path, header = FALSE , stringsAsFactors = FALSE)
  gr <- GRanges(seqnames = Rle(bed_data[[1]]),
                ranges = IRanges(start = bed_data[[2]] + 1, end = bed_data[[3]]),
                strand = ifelse(length(bed_data) >= 6, bed_data[[6]], "*"))
  return(gr)
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
summary <- data.frame()
for(i in c(1:length(tissues))){
  tissue_a <- tissues[i]
  gr1 <- read_bed_to_granges(paste0("data/samples/",tissue_a,"/H3K27me3/bed/H3K27me3_young_old_merge-W5000-G10000-E100.bed"))
  for(j in c(1:length(tissues))){
    tissue_b <- tissues[j]
    gr2 <- read_bed_to_granges(paste0("data/samples/",tissue_b,"/H3K27me3/bed/H3K27me3_young_old_merge-W5000-G10000-E100.bed"))
    combined_gr <- c(gr1, gr2)
    reduced_gr <- as.data.frame(reduce(combined_gr))
    overlaps <- findOverlaps(gr1, gr2)
    overlap_gr <- as.data.frame(reduce(c(gr1[queryHits(overlaps)], gr2[subjectHits(overlaps)])))
    t_summary <- data.frame(tissue1=tissue_label_change(tissue_a),tissue2=tissue_label_change(tissue_b),
                            overlap=nrow(overlap_gr), union=nrow(reduced_gr),percent=nrow(overlap_gr)/nrow(reduced_gr)*100)
    summary <- rbind(summary,t_summary)
    }
}
to_plot <- reshape2::dcast(summary,tissue1~tissue2,value.var="percent")
rownames(to_plot) <- to_plot$tissue1
to_plot <- to_plot[,-1]
breaks <- c(seq(50, 70, length.out = 40), seq(71, 85, length.out = 20), seq(86, 100, length.out = 40))
pheatmap::pheatmap(to_plot,breaks = breaks)
