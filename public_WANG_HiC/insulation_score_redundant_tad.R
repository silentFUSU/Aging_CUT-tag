rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)

resolution <- "20000"
samples <- c("G1","G1","DS1","DS2")
boundary_list <- data.frame()
for(sample in samples){
  boundary <- read.table(paste0("data/public_data/WANG_cellular_aging_HiC/TAD/insulation_score/",sample,"/",sample,"_",resolution,"_dense.is500001.ids200001.insulation.boundaries"))
  boundary <- boundary[,c(7:8)]
  boundary <- boundary %>%  
    separate(V7, into = c("bin", "org", "pos"), sep = "\\|") %>%  
    separate(pos, into = c("chr", "start_end"), sep = ":") %>%  
    separate(start_end, into = c("start", "end"), sep = "-")  
  boundary <- boundary[,-c(1:2)]
  colnames(boundary)[4] <- "boundaryScore"
  boundary$sample <- sample
  boundary$chr <- factor(boundary$chr, levels = paste0("chr",c(1:22,"X","Y")))
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
  while(sorted_boundary_list[i,"chr"]==sorted_boundary_list[j,"chr"] &  sorted_boundary_list[j,"start"]-sorted_boundary_list[i,"start"] <= 60000 & j <=nrow(boundary_list) ){
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
write.csv(redundant_boundary,"data/public_data/WANG_cellular_aging_HiC/TAD/insulation_score/redundant_TAD.csv",row.names = F)
redundant_boundary$label <- paste0(redundant_boundary$chr,":",redundant_boundary$start,"-",redundant_boundary$end)
redundant_boundary_list <- paste0(redundant_boundary$chr,":",redundant_boundary$start,"-",redundant_boundary$end)
bedpe <- data.frame()
for(row_num in c(1:(nrow(redundant_boundary)-1))){
  t_bedpe <- data.frame(chr=redundant_boundary[row_num,1],
                      start = (redundant_boundary[row_num,2]+redundant_boundary[row_num,3])/2, 
                      end = (redundant_boundary[row_num+1,2]+redundant_boundary[row_num+1,3])/2)
  bedpe <- rbind(bedpe,t_bedpe)  
}
    

bedpe <- data.frame(chr1=bedpe$chr, x1=bedpe$start, x2=bedpe$end, 
                    chr2=bedpe$chr, y1=bedpe$start, y2=bedpe$end)
write.table(bedpe,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/all_samples_",resolution,"_redundant_tads.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)

