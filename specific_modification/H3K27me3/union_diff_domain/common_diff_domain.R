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

tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
antibody <- "H3K27me3"
common_diff_peak <- function(antibody,tissues){
  domain <- data.frame(Geneid = character(),  
                   Chr = character(),
                   Start = numeric(),
                   End = numeric(),
                   Significant = character(),  
                   tissue = character(),
                   stringsAsFactors = FALSE) 
  for (i in c(1:length(tissues))){
    tissue <- tissues[i]
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_edd_domain_merged_diff_after_remove_batch_effect.csv")) 
    df <- df[which(df$Significant!="Stable"),c("Geneid","Chr","Start","End","Significant")]
    if(nrow(df) >0){
      df$tissue <- tissue
      df <- unique(df)
      domain <- rbind(domain,df)
    }
  }
  domain_file <- read.table(paste0("data/samples/all/",antibody,"/bed/",antibody,"_edd_domain_merged.bed"))
  domain_file$Geneid <- paste0(domain_file$V1,":",domain_file$V2,"-",domain_file$V3)
  colnames(domain_file)[4]<-"Geneid"
  
  increase <- domain[which(domain$Significant=="Up"),]
  increase_count <- increase %>%   
    count(Geneid)
  increase_tissue <- increase %>%   
    group_by(Geneid) %>%   
    summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  increase_count <- merge(increase_count,increase_tissue,by="Geneid")
  increase_count <- merge(increase_count,domain_file,by="Geneid")
  colnames(increase_count)[4:6] <- c("chr","start","end")
  
  decrease <- domain[which(domain$Significant=="Down"),]
  decrease_count <- decrease %>%   
    count(Geneid)
  decrease_tissue <- decrease %>%   
    group_by(Geneid) %>%   
    summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  decrease_count <- merge(decrease_count,decrease_tissue,by="Geneid")
  decrease_count <- merge(decrease_count,domain_file,by="Geneid")
  colnames(decrease_count)[4:6] <- c("chr","start","end")
  write.csv(increase_count,paste0("data/samples/all/",antibody,"/common_increase_edd_domain_merged_after_remove_batch_effect.csv"),row.names = F)
  write.csv(decrease_count,paste0("data/samples/all/",antibody,"/common_decrease_edd_domain_merged_after_remove_batch_effect.csv"),row.names = F)
  
}
