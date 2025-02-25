rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
tissue <- "thymus"
resolution <- "10000"
tad_format_transfer <- function(tissue,resolution){
  chromosomes <- c(paste0("chr",c(1:19,"X","Y")))
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    chr_list <- list()
    for(chr in chromosomes){
      df <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"/",sample,"_",resolution,"_",chr,"_dense.is500001.ids200001.insulation.boundaries.bed"),skip=1)
      chr_df <- data.frame(chr=as.character(),start=as.numeric(),end=as.numeric())
      for(j in c(1:(nrow(df)-1))){
        t_chr_df <- data.frame(chr=chr,
                               start = (df[j,2]+df[j,3])/2, 
                               end = (df[j+1,2]+df[j+1,3])/2)
        chr_df <- rbind(chr_df,t_chr_df)
      }
      chr_list[[chr]]<-chr_df
    }
    tads <- Reduce(function(x, y) rbind(x, y),chr_list)
    write.csv(tads,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"_",resolution,"_tads.csv"),row.names = F)
  }
}
# tissues <- c("colon","kidney","liver","lung","CB","brain")
tissues <- c("bonemarrow","stomach","heart")
for(tissue in tissues){
  tad_format_transfer(tissue,resolution)
}

convert2bedpe <- function(tissue, resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(sample in search_table$sample_name){
    df <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"_",resolution,"_tads.csv"))
    bedpe <- data.frame(chr1=df$chr,x1=df$start,x2=df$end,
                        chr2=df$chr,y1=df$start,y2=df$end)
    write.table(bedpe,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"_",resolution,"_tads.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
    }
}
for(tissue in tissues){
  convert2bedpe(tissue,resolution)
}
