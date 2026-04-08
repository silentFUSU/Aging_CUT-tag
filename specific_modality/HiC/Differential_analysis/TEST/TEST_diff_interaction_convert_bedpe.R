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
resolution <- "20000"
convert2bedpe <- function(tissue,resolution){
  df <-  read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff_larger_250000.csv"))
  df$region1 <- df$X
  df$region2 <- df$X
  df <- df[order(df$FDR.old.young),]
  df <- df[,c("region1","region2","Significant")]
  df <- df %>%
    separate(region1, into = c("chr1", "x1", "x2"), sep = "-", convert = TRUE)
  df <- df %>%
    separate(region2, into = c("chr2", "y1", "y2"), sep = "-", convert = TRUE)
  
  write.table(head(df[which(df$Significant == "Up"), c(1:6)], 10),paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff_larger_250000_increase.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names =T)
  write.table(head(df[which(df$Significant == "Down"), c(1:6)], 10),paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff_larger_250000_decrease.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)
}
tissues <- sort(c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus","skin","muscle","cecum","ileum","spleen","pancreas"))
for(tissue in tissues){
  convert2bedpe(tissue,resolution)
}

