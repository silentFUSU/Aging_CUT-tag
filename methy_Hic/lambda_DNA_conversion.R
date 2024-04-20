rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)  
file_path <- args[1]
# file_path <- "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/methy_HiC/20240411_hundep2_methylHic_data/WJH_MethylHic_10/vcf/WJH_MethylHic_10.chrL.calmd.cpg.raw.sort.vcf"
data <- read.table(file_path, header=F, comment.char="#") 
data <- data %>%  
  separate(V8, into = paste0("V8_", 1:5), sep = ";", extra = "merge")
data <- data[which(data$V8_1=="CS=+"),]

data <- data %>% 
  separate(V8_3, into = paste0("V8_3_", 1:2), sep = "=", extra = "merge")
data$V8_3_2 <- as.numeric(data$V8_3_2)

data <- data %>%
  separate(V10, into = paste0("V10_", 1:11), sep = ":", extra = "merge")

data_noNA <- data %>% filter(!is.na(V10_4), V10_4 != ".") 
data_noNA$V10_4 <- as.numeric(data_noNA$V10_4)
total <- sum(data_noNA$V8_3_2)
noconvert <- sum(data_noNA$V10_4)
convert <- 1-(noconvert/total)
trimmed_path <- sub("(.*)chrL.*", "\\1", file_path) 
fileConn <- file(paste0(trimmed_path,"calmd.nodup.cpg.raw.vcf.MethySummarizeList.txt"),open = "a" ) 
writeLines(paste0("Lambda DNA conversation rate: ",convert*100,"%"),fileConn)
close(fileConn) 
