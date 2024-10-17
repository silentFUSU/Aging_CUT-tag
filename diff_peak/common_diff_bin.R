rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(tidyr)
library(UpSetR)
library(ChIPseeker)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT")
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
antibody <- "H3K9me3"
commom_diff_bin <- function(antibody,tissues){
  bin<-data.frame(Geneid = character(),  
                    Chr = character(),
                    Start = numeric(),
                    End = numeric(),
                    Significant = character(),  
                    tissue = character(),
                    stringsAsFactors = FALSE) 
  for (i in c(1:length(tissues))){
    tissue <- tissues[i]
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff.csv")) 
    df <- df[which(df$Significant!="Stable"),c("Geneid","Chr","Start","End","Significant")]
    if(nrow(df) >0){
      df$tissue <- tissue
      df <- unique(df)
      bin <- rbind(bin,df)
    }
  }
  bin_file <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_",bin_size(antibody),"_bins.bed"))
  colnames(bin_file)[4]<-"Geneid"
  
  increase <- bin[which(bin$Significant=="Up"),]
  increase_count <- increase %>%   
    count(Geneid)
  increase_tissue <- increase %>%   
    group_by(Geneid) %>%   
    summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  increase_count <- merge(increase_count,increase_tissue,by="Geneid")
  increase_count <- merge(increase_count,bin_file,by="Geneid")
  colnames(increase_count)[4:6] <- c("chr","start","end")
  
  decrease <- bin[which(bin$Significant=="Down"),]
  decrease_count <- decrease %>%   
    count(Geneid)
  decrease_tissue <- decrease %>%   
    group_by(Geneid) %>%   
    summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  decrease_count <- merge(decrease_count,decrease_tissue,by="Geneid")
  decrease_count <- merge(decrease_count,bin_file,by="Geneid")
  write.csv(increase_count,paste0("data/samples/all/",antibody,"/common_increase_",bin_size(antibody),"_bins.csv"))
  write.csv(decrease_count,paste0("data/samples/all/",antibody,"/common_decrease_",bin_size(antibody),"_bins.csv"))
  
  }
 




