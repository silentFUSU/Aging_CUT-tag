rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(patchwork)
library(stringr)
library(dplyr)
library(tidyr)
library(UpSetR)
library(ChIPseeker)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","FC","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum")

genes<-data.frame(X = character(),  
                  Significant = character(),  
                  tissue = character(),
                  stringsAsFactors = FALSE)  

for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_nodup.csv")) 
  df <- df[which(df$Significant!="Stable"),c("X","Significant")]
  if(nrow(df) >0){
    df$tissue <- tissue
    genes <- rbind(genes,df)
  }
}
colnames(genes)[1] <-"Geneid"
increase <- genes[which(genes$Significant=="Up"),]
increase_count <- increase %>%   
  count(Geneid)
increase_tissue <- increase %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
increase_count <- merge(increase_count,increase_tissue,by="Geneid")
increase_count <- merge(increase_count,bin_file,by="Geneid")

colnames(genes)[1] <-"Geneid"
decrease <- genes[which(genes$Significant=="Down"),]
decrease_count <- decrease %>%   
  count(Geneid)
decrease_tissue <- decrease %>%   
  group_by(Geneid) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
decrease_count <- merge(decrease_count,decrease_tissue,by="Geneid")

