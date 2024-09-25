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
library(DSS)

args <- commandArgs(trailingOnly = TRUE)  
if (length(args) < 1) {  
  stop("No tissue argument provided")  
}  
tissue <- args[1]  
print(paste("Tissue is:", tissue))  

bdg2dsstxt <- function(tissue){
  data_path <- paste0("data/samples/WGBS/",tissue,"/")
  dir.create(paste0(data_path,"DSS_table"))
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  for(sample in search_table$sample_name){
    df <- fread(paste0(data_path,"bdg/",sample,"_CpG.bdg"),sep="\t")
    df <- df[,c(1,2,5,4)]  
    colnames(df) <- c("chr", "pos", "N", "X")
    fwrite(df,paste0(data_path,"/DSS_table/",sample,".txt"), sep = "\t")  
    }
}
bdg2dsstxt(tissue)

DSS_DMR <- function(tissue){
  data_path <- paste0("data/samples/WGBS/",tissue,"/")
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  young <- search_table$sample_name[which(search_table$age=="3M")]
  old <- search_table$sample_name[which(search_table$age=="24M")]
  young_list <- list()
  age_young <- vector()
  for(i in c(1:length(young))){
    sample <- young[i]
    young_list[[i]] <- fread(paste0(data_path,"/DSS_table/",sample,".txt"),sep="\t")
    names(young_list)[i] <- sample
    age_young <- append(age_young,paste0(search_table$age[which(search_table$sample_name == sample)],"_rep",i))
  }
  

  old_list <- list()
  age_old <- vector()
  for(i in c(1:length(old))){
    sample <- old[i]
    old_list[[i]] <- fread(paste0(data_path,"/DSS_table/",sample,".txt"),sep="\t")
    names(old_list)[i] <- sample
    age_old <- append(age_old,paste0(search_table$age[which(search_table$sample_name == sample)],"_rep",i))
  }
  
  merge_list <- append(young_list,old_list)
  BSobj = makeBSseqData( merge_list,append(age_young,age_old))
  dmlTest.sm = DMLtest(BSobj, group1=age_old, group2=age_young, smoothing=TRUE,ncores=10)                   
  # dmls <- callDML(dmlTest.sm, p.threshold=.01, delta=0.1)
  dmls <- callDML(dmlTest.sm, p.threshold=.01, delta=0)
  write.table(dmls, paste0(data_path,"/DSS_table/",tissue,"_DML_delta0.txt"), row.names=F, sep='\t', quote=F)
  
  # dmrs <- callDMR(dmlTest.sm, p.threshold=.01, delta=0.1)
  dmrs <- callDMR(dmlTest.sm, p.threshold=.01, delta=0)
  write.table(dmrs, paste0(data_path,"/DSS_table/",tissue,"_DMR_delta0.txt"), row.names=F, sep='\t', quote=F)
  
  saveRDS(BSobj,paste0(data_path,"/DSS_table/",tissue,"_Bsobj_delta0.rds"))
  saveRDS(dmlTest.sm,paste0(data_path,"/DSS_table/",tissue,"_dmlTest_delta0.rds"))
  # increase <- DMR[which(DMR$areaStat > 0),c("chr","start","end")]
  # write.table(increase, file=paste0(data_path,"/DSS_table/bed/DMR_increase.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  # decrease <- DMR[which(DMR$areaStat < 0),c("chr","start","end")]
  # write.table(decrease, file=paste0(data_path,"DSS_table/bed/DMR_decrease.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
  }
DSS_DMR(tissue)
