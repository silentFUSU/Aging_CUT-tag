rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)  
library(parallel)  
tss_ref <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
tss_ref$V2 <- tss_ref$V2 - 1000
tss_ref$V3 <- tss_ref$V3 + 1000
tss_ref <- tss_ref[which(tss_ref$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
tss_ref <- tss_ref[,c(1:3,6)]
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")


compress_to_TSS_1kb <- function(tissue){
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
    t_tss_ref <- data.table(tss_ref)
    
    setDT(df)  
    setDT(t_tss_ref)  
    setkey(df, V1, V2, V3) 
    setkey(t_tss_ref, V1, V2, V3) 
    overlaps <- foverlaps(df, t_tss_ref, type = "any", nomatch = 0L) 
    results <- overlaps[, .(total_V4 = sum(V4), total_V5 = sum(V5)), by = .(V6)]  
    results$percent <- results$total_V4/results$total_V5
    colnames(results)[1] <- "gene_name"
    results$tissue <- tissue
    results$sample <- sample
    bin_summary <- rbind(bin_summary,results)
  }
  write.csv(bin_summary,paste0("data/samples/WGBS/",tissue,"/compress2bin/TSS_region_1kb_all_depth.csv"),row.names = F)
}
tissues <- c("liver","lung","mammarygland","kidney","Hip","ileum")
for(tissue in tissues){
  compress_to_TSS_1kb(tissue)
}

