rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(dplyr)
library(reshape2)
library(VennDiagram)
df <- read.csv("data/samples/muscle/H3K36me3/H3K36me3_10kb_bins_diff_after_remove_batch_effect.csv")
ggplot() +
  geom_point(data=df[which(df$Significant=="Stable"),], mapping=aes(logCPM,LogFC.old.young), size=2, color="grey") +
  geom_point(data=df[which(df$Significant=="Up"),], mapping=aes(logCPM, LogFC.old.young), size=2, color="red") +
  geom_point(data=df[which(df$Significant=="Down"),], mapping=aes(logCPM,LogFC.old.young), size=2, color="blue") +
labs(x="log2(CPM)",
     y="log2(FC)") +
  theme(legend.position = "bottom",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme_bw()
H3K36me3_yong_peak <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_young_merge-W1000-G3000-E100.bed")
H3K36me3_yong_peak_down <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_10kb_bins_diff_after_remove_batch_effect_down_in_young_peaks_unique_sort.bed")
H3K36me3_yong_peak_up <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_10kb_bins_diff_after_remove_batch_effect_up_in_young_peaks_unique_sort.bed")
dir.create("result/muscle/H3K36me3/")
venn.diagram(
  x = list(H3K36me3_yong_peak$V4,H3K36me3_yong_peak_down$V4,H3K36me3_yong_peak_up$V4),
  category.names = c("Young samples all peaks" , "Decrease peaks(old-young) " , "Increas peaks(old-young)"),
  filename = 'result/muscle/H3K36me3/bins_up_down_in_young_peaks_venn.png',
  output=TRUE
)

H3K36me3_old_peak <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_old_merge-W1000-G3000-E100.bed")
H3K36me3_old_peak_down <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_10kb_bins_diff_after_remove_batch_effect_down_in_old_peaks_unique_sort.bed")
H3K36me3_old_peak_up <- read.table("data/samples/muscle/H3K36me3/bed/H3K36me3_10kb_bins_diff_after_remove_batch_effect_up_in_old_peaks_unique_sort.bed")
dir.create("result/muscle/H3K36me3/")
venn.diagram(
  x = list(H3K36me3_old_peak$V4,H3K36me3_old_peak_down$V4,H3K36me3_old_peak_up$V4),
  category.names = c("old samples all peaks" , "Decrease peaks(old-young) " , "Increas peaks(old-young)"),
  filename = 'result/muscle/H3K36me3/bins_up_down_in_old_peaks_venn.png',
  output=TRUE
)
