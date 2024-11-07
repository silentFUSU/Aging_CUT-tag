rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
library(tidyr)
tissue <- "lung"
resolution <- "10000"
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

loop_combine <- function(tissue, resolution){
  chromosomes <- paste0("chr",c(1:19,"X","Y"))
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    chr_list <- list()
    for(chr in chromosomes){
      if(file.exists(paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/",sample,"/",chr,"/merged_loops.bedpe")) &&  file.info(paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/",sample,"/",chr,"/merged_loops.bedpe"))$size > 0){
        loop <- read.table(paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/",sample,"/",chr,"/merged_loops.bedpe"))
        loop <- loop[,c(1:6)]
        chr_list[[chr]] <- loop
      }
    }
    loop_summary <- Reduce(function(x, y) rbind(x,y), chr_list) 
    loop_summary$label <- paste0(loop_summary$V1,"-",loop_summary$V2,"-",loop_summary$V3,"_",loop_summary$V4,"-",loop_summary$V5,"-",loop_summary$V6)
    write.csv(loop_summary,paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/",sample,"/",sample,"_",resolution,"_loop.csv"),row.names = F)
  }
}
loop_combine(tissue,resolution)

loop_change <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  sample_list <- list(young=search_table$sample_name[which(search_table$age=="3M")],
                      old=search_table$sample_name[which(search_table$age=="24M")])
  ages <- c("young","old")
  loop_summary <- data.frame()
  loop_list <- list(young=list(), old=list())
  for(age in ages){
    for(sample in sample_list[[age]]){
      df <- read.csv(paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/",sample,"/",sample,"_",resolution,"_loop.csv"))
      t_loop_summary <- data.frame(sample=sample,loops=nrow(df),age=age)  
      loop_summary <- rbind(loop_summary,t_loop_summary)
      loop_list[[age]][[sample]] <- df
    }
  }
  young_list <- lapply(loop_list[["young"]], function(df) df$label)
  young_loop <- Reduce(intersect, young_list)
  
  old_list <- lapply(loop_list[["old"]], function(df) df$label)
  old_loop <- Reduce(intersect, old_list)
}
