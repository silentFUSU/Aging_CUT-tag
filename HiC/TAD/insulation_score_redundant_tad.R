rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)
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
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "cecum"
resolution <- "20000"
insulation_redundant_TAD <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  boundary_list <- data.frame()
  for(sample in search_table$sample_name){
    boundary <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"/",sample,"_",resolution,"_dense.is500001.ids200001.insulation.boundaries"))
    boundary <- boundary[,c(7:8)]
    boundary <- boundary %>%  
      separate(V7, into = c("bin", "org", "pos"), sep = "\\|") %>%  
      separate(pos, into = c("chr", "start_end"), sep = ":") %>%  
      separate(start_end, into = c("start", "end"), sep = "-")  
    boundary <- boundary[,-c(1:2)]
    colnames(boundary)[4] <- "boundaryScore"
    boundary$sample <- sample
    boundary$chr <- factor(boundary$chr, levels = paste0("chr",c(1:19,"X","Y")))
    boundary$start <- as.numeric(boundary$start)
    boundary$end <- as.numeric(boundary$end)
    boundary_list <- rbind(boundary_list,boundary)
  }
  sorted_boundary_list <- boundary_list %>%
    arrange(chr, start, boundaryScore)  
  sorted_boundary_list$choose_or_not <- "No"
  i=1
  redundant_boundary <- data.frame()
  while(i < nrow(boundary_list)){
    j=i+1
    while(sorted_boundary_list[i,"chr"]==sorted_boundary_list[j,"chr"] &  sorted_boundary_list[j,"start"]-sorted_boundary_list[i,"start"] <= 50000 & j <=nrow(boundary_list) ){
        if(sorted_boundary_list[i,"boundaryScore"] < sorted_boundary_list[j,"boundaryScore"]){
          i = j
          j+1
        }else{
          j=j+1
        }
    }
    sorted_boundary_list[i,"choose_or_not"] <- "Yes"
    redundant_boundary <- rbind(redundant_boundary,sorted_boundary_list[i,])
    i=j
  }
  
  redundant_boundary <- redundant_boundary[,-ncol(redundant_boundary)]
  write.csv(redundant_boundary,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD.csv"),row.names = F)
  redundant_boundary$label <- paste0(redundant_boundary$chr,":",redundant_boundary$start,"-",redundant_boundary$end)
  redundant_boundary_list <- paste0(redundant_boundary$chr,":",redundant_boundary$start,"-",redundant_boundary$end)

  bedpe <- data.frame(chr=as.character(),start=as.numeric(),end=as.numeric())
  i=2
  while(i <= nrow(redundant_boundary)){
    if(redundant_boundary[i,"chr"] == redundant_boundary[(i-1),"chr"]){
      t_bedpe <- data.frame(chr=redundant_boundary$chr[i],start=redundant_boundary$start[i-1],end=redundant_boundary$start[i]) 
      bedpe <- rbind(bedpe,t_bedpe)
      i <- i+1
    }else{
      i <- i+1
    }
  }

  bedpe <- data.frame(chr1=bedpe$chr, x1=bedpe$start, x2=bedpe$end, 
                      chr2=bedpe$chr, y1=bedpe$start, y2=bedpe$end)
  write.table(bedpe,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/all_samples_",resolution,"_redundant_tads.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  }
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
#bonemarrow heart stomach mammarygland
for(tissue in tissues){
  insulation_redundant_TAD(tissue,resolution)
}


