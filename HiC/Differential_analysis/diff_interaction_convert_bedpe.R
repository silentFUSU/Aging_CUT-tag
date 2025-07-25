rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
tissue <- "lung"
resolution <- "200000"
convert2bedpe <- function(tissue,resolution){
  re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
  re$Significant <- "Stable"
  re$Significant[which((re$V5 < 0.05 & re$V6 < 0.05) & re$V4 < 0)] <- "Down"
  re$Significant[which((re$V5 < 0.05 & re$V6 < 0.05) & re$V4 > 0)] <- "Up"
  re_sig <- re[which(re$Significant!="Stable"),]
  
  re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4 & abs(re_sig$V2 - re_sig$V3) <= 200),]
  re_sig$V1 <- paste0("chr",re_sig$V1)
  re_sig$V1[which(re_sig$V1=="chr20")] <- "chrX"
  re_sig$label1 <- paste(re_sig$V1,re_sig$V2,sep = "-")
  re_sig$label2 <- paste(re_sig$V1,re_sig$V3,sep = "-")

  
  HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed <- bed[which(bed$V1 %in% paste0("chr",c(1:19,"X"))),]
  bed$V5 <- "NA"
  count <- 1
  bed[1,"V5"] <- count
  for(i in c(2:nrow(bed))){
    if(bed[i,"V1"] != bed[i-1,"V1"]){
      count <- 1
      bed[i,"V5"] <- count
    }else{
      count <- count+1
      bed[i,"V5"] <- count
    }
  }
  bed$region <- paste(bed$V1,bed$V2,bed$V3,sep = "-")
  bed$label <- paste(bed$V1,bed$V5,sep = "-")
  
  re_sig <- merge(re_sig,bed,by.x="label1",by.y="label")
  colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region1"
  re_sig <- merge(re_sig,bed,by.x="label2",by.y="label")
  colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region2"
  
  df <- re_sig[,c("region1","region2","Significant")]
  df <- df %>%
    separate(region1, into = c("chr1", "x1", "x2"), sep = "-", convert = TRUE)
  df <- df %>%
    separate(region2, into = c("chr2", "y1", "y2"), sep = "-", convert = TRUE)
  write.table(df[which(df$Significant=="Up"),c(1:6)],paste0("data/samples/HiC/",tissue,"/differential_analysis/Wang_output_",resolution,"_increase.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(df[which(df$Significant=="Down"),c(1:6)],paste0("data/samples/HiC/",tissue,"/differential_analysis/Wang_output_",resolution,"_decrease.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
for(tissue in tissues){
  convert2bedpe(tissue,resolution)
}

