rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(enrichplot)

options(scipen =0)  
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
genes <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
genes <- genes[,c("V1","V2","V3","V6")]
genes <- as.data.table(genes)
setDT(genes)
setkey(genes,V1,V2,V3)
kmean <- "kmeans1"
regions <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_chr <- strsplit(as.character(regions$label), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
regions <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=regions$cluster)
regions$start <- as.numeric(regions$start)
regions$end <- as.numeric(regions$end)
regions <-regions[which(regions$cluster==kmean),]
regions <- as.data.table(regions)
setDT(regions)
setkey(regions,chr,start,end)

overlaps <-  foverlaps(genes, regions, type = "any", nomatch = 0L)  
genes <- unique(overlaps$V6)
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_summary <- data.frame()
genes_up <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  df <- df[which(df$X %in% genes),]
  t_tissue_summary <- as.data.frame(table(df$Significant))
  t_tissue_summary <- rbind(t_tissue_summary,data.frame(Var1="Filtered out", Freq=length(genes)-sum(t_tissue_summary$Freq)))
  t_tissue_summary$percent <- t_tissue_summary$Freq / sum(t_tissue_summary$Freq)
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
  t_genes_up <- df[which(df$Significant=="Up"),"X",drop=F]
  t_genes_up$tissue <- tissue_label_change(tissue)
  genes_up <- rbind(genes_up,t_genes_up)
}

genes_up_count <- genes_up %>%   
  count(X)
genes_up_tissue <- genes_up %>%   
  group_by(X) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
genes_up_count <- merge(genes_up_count,genes_up_tissue,by="X")


