rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
tissue <- "BAT"
search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
H3K27me3_change_in_H3K9me3_decrease_region <- function(tissue){
  H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K9me3_decrease_region <- H3K9me3[which(H3K9me3$Significant_bar=="Down"),]
  write.table(H3K9me3_decrease_region[,c("Chr","Start","End","Geneid")], 
              paste0("data/samples/BAT/H3K9me3/bed/H3K9me3_10kb_bins_diff_after_remove_batch_effect_down.bed"), 
              quote = F, row.names = F, col.names = F)
  search_table <- search_table[which(search_table$tissue == tissue & search_table$antibody== "H3K27me3"),]
  young <-search_table$sample_name[which(search_table$age=="3m")]
  old <- search_table$sample_name[which(search_table$age=="24m")]
  dir.create(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/matrix"))
  dir.create(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/plot"))
  }