rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999)  
DMR_split <- function(tissue){
  DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01.txt"),header = T)
  dir.create(paste0("data/samples/WGBS/",tissue,"/DSS_table/bed/"))
  increase <- DMR[which(DMR$areaStat > 0),c("chr","start","end")]
  # increase <- increase[order(increase$areaStat,decreasing = T),]
  # increase <- increase[c(1:50000),c("chr","start","end")]
  write.table(increase, file=paste0("data/samples/WGBS/",tissue,"/DSS_table/bed/",tissue,"_DMR_increase_delta01.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  decrease <- DMR[which(DMR$areaStat < 0),c("chr","start","end")]
  write.table(decrease, file=paste0("data/samples/WGBS/",tissue,"/DSS_table/bed/",tissue,"_DMR_decrease_delta01.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver","lung","muscle","pancreas","skin","spleen","stomach","testis","thymus","tongue","iWAT","mammarygland","ovary","uterus")
for( tissue in tissues) {
  DMR_split(tissue)
}

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver","lung","muscle","skin","spleen","stomach","testis","thymus","tongue","iWAT","mammarygland","uterus")
for( tissue in tissues) {
  t_DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01.txt"),header = T)
  t_DMR$label <- paste0(t_DMR$chr,t_DMR$start,t_DMR$end,sep="-")
  df <- df[which(!df %in% t_DMR$label)]
  }
