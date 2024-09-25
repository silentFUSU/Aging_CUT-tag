rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)  
library(parallel)  
gene_body_ref <- read.table("~/ref_data/for_normal_mapping/TSS/refBed/mm10_gene.bed")
gene_body_ref <- gene_body_ref[which(gene_body_ref$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
gene_body_ref <- gene_body_ref[,c(1:4)]
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")


compress_to_gene_body <- function(tissue){
  samples <- search_table$sample_name[which(search_table$tissue==tissue)]
  bin_summary <- data.frame(gene_name = character(),
                            total_V4 = numeric(),
                            total_V5 = numeric(),
                            percent = numeric(),
                            tissue = character(),
                            sample = character())
  for(sample in samples){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    df$V1 <- factor(df$V1,paste0("chr",c(1:19,"X","Y")))
    df$V3 <- df$V2
    df <- df[order(df$V1,df$V2),]
    t_gene_body_ref <- data.table(gene_body_ref)
    
    setDT(df)  
    setDT(t_gene_body_ref)  
    setkey(df, V1, V2, V3) 
    setkey(t_gene_body_ref, V1, V2, V3) 
    overlaps <- foverlaps(df,t_gene_body_ref, type = "any", nomatch = 0L) 
    results <- overlaps[, .(total_V4 = sum(i.V4), total_V5 = sum(V5)), by = .(V4)]  
    results$percent <- results$total_V4/results$total_V5
    colnames(results)[1] <- "gene_name"
    results$tissue <- tissue
    results$sample <- sample
    bin_summary <- rbind(bin_summary,results)
  }
  write.csv(bin_summary,paste0("data/samples/WGBS/",tissue,"/compress2bin/gene_body_region_all_depth.csv"),row.names = F)
}
tissues <- c("liver","lung","mammarygland","kidney","Hip","ileum")
for(tissue in tissues){
  compress_to_gene_body(tissue)
}

