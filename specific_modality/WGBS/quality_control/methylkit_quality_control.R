rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(methylKit)
library(data.table)
tissue <- "lung"
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
methykit_data <- function(tissue){
  data_path <- "data/samples/WGBS/"
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  dir.create(paste0(data_path,"/",tissue,"/methykit_table/"))
  for(i in c(1:length(t_search_table$sample_name))){
    df <- fread(paste0(data_path,tissue,"/bdg/",t_search_table$sample_name[i],"_CpG.bdg"),sep = "\t")
    df$chrBase <- paste0(df$V1,".",df$V2)
    df$V6[which(df$V6=="+")] <- "F"
    df$V6[which(df$V6=="-")] <- "R"
    methykit_df <- data.frame(chrBase=df$chrBase,
                              chr=df$V1,
                              base=df$V2,
                              strand=df$V6,
                              coverage=df$V5,
                              freqC=round(df$V4/df$V5,4)*100,
                              freqT=100-round(df$V4/df$V5,4)*100)
    fwrite(methykit_df, file = paste0(data_path, "/", tissue, "/methykit_table/", t_search_table$sample_name[i], "_CpG.txt"))
  }
}
methykit_data(tissue)
methykit_QC <- function(tissue){
  methyl_obj <- methRead(data_list,
                         sample.id=list("test1","test2","ctrl1","ctrl2"),
                         assembly="hg18",
                         treatment=c(1,1,0,0),
                         context="CpG",
                         dbtype = "tabix",
                         dbdir = "methylDB"
  )
}