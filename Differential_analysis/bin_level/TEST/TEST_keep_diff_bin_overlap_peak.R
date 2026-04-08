rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff.csv"))
    bin_in_peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_young_old_merge-W1000-G3000-E100.bed"))
    diff_bin_in_peak <- diff[which(diff$Geneid %in% bin_in_peak$V4),]
    write.csv(diff_bin_in_peak,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_overlap_young_old_peaks.csv"))
  }
}

antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff.csv"))
    bin_in_peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_1kb_in_young_old_merge_macs_narrowpeak.bed"))
    diff_bin_in_peak <- diff[which(diff$Geneid %in% bin_in_peak$V4),]
    write.csv(diff_bin_in_peak,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff_overlap_young_old_peaks.csv"))
  }
}

antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody<-antibodys[j]
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
    bin_in_peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_macs_narrowpeak.bed"))
    diff_bin_in_peak <- diff[which(diff$Geneid %in% bin_in_peak$V4),]
    write.csv(diff_bin_in_peak,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins_overlap_peaks_diff_after_remove_batch_effect.csv"))
  }
}
