rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak_num_summary <- data.frame()
peak_annotation_summary <- data.frame()
for(tissue in tissues){
  # df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
  df <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_macs_young_old_narrowpeak_summits_spm3.bed"))
  gr <- GRanges(seqnames = df$V1,
                ranges = IRanges(start = df$V2, end = df$V3))
  peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
  t_peak_num_summary <- data.frame(tissue=tissue_label_change(tissue),num=peakAnno@peakNum)
  peak_num_summary <- rbind(peak_num_summary,t_peak_num_summary)
  t_peak_annotation_summary <- peakAnno@annoStat
  colnames(t_peak_annotation_summary)[2] <- tissue_label_change(tissue)
  if(nrow(peak_annotation_summary)==0){
    peak_annotation_summary <- t_peak_annotation_summary
  }else{
    peak_annotation_summary <- merge(peak_annotation_summary,t_peak_annotation_summary,by="Feature")
  }
}
annotation_long <- peak_annotation_summary %>%
  pivot_longer(-Feature, names_to = "tissue", values_to = "percent")
annotation_long <- annotation_long %>%
  left_join(peak_num_summary, by = "tissue") %>%
  mutate(count = percent / 100 * num)
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,sort(unique(annotation_long$Feature)))
p <- ggplot(annotation_long, aes(x = tissue, y = count, fill = Feature)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values=color) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Feature Count", title = "Peak Annotation Summary by Tissue")
ggsave("result/Sup_figures/ATAC_peak_number_annotation.pdf",p,width = 10,height = 8)
