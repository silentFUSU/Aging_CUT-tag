rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
library(tidyr)
tissue <- "thymus"
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
  resolution_label <- paste0(as.numeric(resolution)/1000,"kb")
  chromosomes <- paste0("chr",c(1:19,"X","Y"))
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    chr_list <- list()
    for(chr in chromosomes){
      if(file.exists(paste0("data/samples/HiC/",tissue,"/loop/mustache/",sample,"_",resolution_label,"/",chr,"/",sample,"_",resolution_label,"_",chr,"_mustache_loop.tsv"))){
        if(length(readLines(paste0("data/samples/HiC/",tissue,"/loop/mustache/",sample,"_",resolution_label,"/",chr,"/",sample,"_",resolution_label,"_",chr,"_mustache_loop.tsv"), warn = FALSE)) > 1){
          loop <- read.table(paste0("data/samples/HiC/",tissue,"/loop/mustache/",sample,"_",resolution_label,"/",chr,"/",sample,"_",resolution_label,"_",chr,"_mustache_loop.tsv"),skip=1)
          loop <- loop[,c(1:6)]
          chr_list[[chr]] <- loop
        }
      }
    }
    loop_summary <- Reduce(function(x, y) rbind(x,y), chr_list) 
    loop_summary$label <- paste0(loop_summary$V1,"-",loop_summary$V2,"-",loop_summary$V3,"_",loop_summary$V4,"-",loop_summary$V5,"-",loop_summary$V6)
    write.csv(loop_summary,paste0("data/samples/HiC/",tissue,"/loop/mustache/",sample,"_",resolution_label,"/",sample,"_",resolution_label,"_mustache_loop.csv"),row.names = F)
  }
}
tissues <- c("bonemarrow","stomach","heart")
for(tissue in tissues){
  loop_combine(tissue,resolution)
}


loop_change <- function(tissue,resolution){
  resolution_label <- paste0(as.numeric(resolution)/1000,"kb")
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  samples_list <- list(young=search_table$sample_name[which(search_table$age=="3M")],
                       old=search_table$sample_name[which(search_table$age=="24M")])
  ages <- c("young","old")
  loop_summary <- data.frame()
  
  for(age in ages){
    for(i in c(1:length(samples_list[[age]]))){
      sample <- samples_list[[age]][i]
      df <- read.csv(paste0("data/samples/HiC/",tissue,"/loop/mustache/",sample,"_",resolution_label,"/",sample,"_",resolution_label,"_mustache_loop.csv"))
      if(nrow(loop_summary)==0){
        loop_summary <- data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df))
      }else{
        loop_summary <- rbind(loop_summary,data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df)))
      }
    }
  }
  loop_summary$age <- factor(loop_summary$age,levels=c("young","old"))
  # color <- read.table("data/samples/7_distinct_color.txt")
  # color <- setNames(color$V1[c(1,3)],c("young","old"))
  p <- ggplot(loop_summary,aes(x=age,y=counts,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text_repel(aes(label = sample), position = position_jitter(width = 0.2),  size = 5) +
    # scale_color_manual(values=color)+
    ggtitle(paste0(tissue_label_change(tissue)," loop counts"))+
    ylim(0,max(loop_summary$counts+1000))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("loop counts")
  print(p)
}
for(tissue in tissues){
  loop_change(tissue,resolution)
}
