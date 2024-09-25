rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
tissues <- c("liver","lung","mammarygland","kidney")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  DML <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DML.txt"),header=T)
  DML <- DML[which(DML$chr %in% c(paste0("chr",c(1:22,"X","Y")))),]
  DML <- DML[which(DML$fdr < 0.05),]
  DML$label <- paste0(DML$chr,"-",DML$pos)
  DML <- DML[,c("label","mu2","mu1")]
  colnames(DML) <- c("label",paste0(tissue,"_young"),paste0(tissue,"_old"))
  if(i == 1){
    DML_summary <- DML
  }else{
    DML_summary <- merge(DML_summary,DML,by="label")
  }
}

rownames(DML_summary) <- DML_summary$label
DML_summary <- DML_summary[,-1]
DML_summary <- DML_summary *100
pheatmap::pheatmap(DML_summary,cluster_cols = F)
