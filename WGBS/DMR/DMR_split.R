rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
DMR_split <- function(tissue){
  DMR <- read.table(paste0("data/raw_data/20240905_DYQ_005-018_WGBS/DSS_table/",tissue,"_DMR.txt"),header = T)
  dir.create(paste0("data/raw_data/20240905_DYQ_005-018_WGBS/DSS_table/bed/"))
  
  increase <- DMR[which(DMR$areaStat > 0),c("chr","start","end")]
  write.table(increase, file=paste0("data/raw_data/20240905_DYQ_005-018_WGBS/DSS_table/bed/",tissue,"_DMR_increase.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  decrease <- DMR[which(DMR$areaStat < 0),c("chr","start","end")]
  write.table(decrease, file=paste0("data/raw_data/20240905_DYQ_005-018_WGBS/DSS_table/bed/",tissue,"_DMR_decrease.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
}
tissue <-"liver"
DMR_split(tissue)
